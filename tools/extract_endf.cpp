// Extracts the fast-neutron data atomsim needs from an ENDF-6 evaluation.
//
//   extract_endf <evaluation.endf> <output.xs>
//
// Writes, for incident energies 25 keV - 20 MeV:
//   XS  MF3 cross-sections (elastic, inelastic levels and continuum, (n,2n),
//       fission, capture) with their Q-values
//   NU  MF1 MT452 total nu-bar(E)
//   ANG MF4 MT2 elastic angular distributions in the centre-of-mass frame,
//       as 32 equiprobable bins of mu per incident energy
// Source: ENDF/B-VIII.0 (Brown et al., Nucl. Data Sheets 148, 1, 2018).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kMinEnergy = 2.5e4, kMaxEnergy = 2.0e7;  // eV
constexpr int kBins = 32;

using Points = std::vector<std::pair<double, double>>;
using Lines = std::vector<std::string>;

// ENDF floats omit the 'e': "1.234567+5" or "-2.5-3".
double number(const std::string& field) {
  std::size_t a = field.find_first_not_of(' '), b = field.find_last_not_of(' ');
  if (a == std::string::npos) return 0.0;
  std::string s = field.substr(a, b - a + 1);
  for (std::size_t i = s.size() - 1; i > 0; --i)
    if ((s[i] == '+' || s[i] == '-') && s[i - 1] != 'e' && s[i - 1] != 'E') {
      s.insert(i, "e");
      break;
    }
  return std::strtod(s.c_str(), nullptr);
}

// Lines of each (MF, MT) section, in file order.
std::map<std::pair<int, int>, Lines> sections(const char* path) {
  std::ifstream in(path);
  if (!in) {
    std::fprintf(stderr, "extract_endf: cannot open %s\n", path);
    std::exit(1);
  }
  std::map<std::pair<int, int>, Lines> out;
  std::string line;
  while (std::getline(in, line)) {
    if (line.size() < 75) continue;
    int mf = std::atoi(line.substr(70, 2).c_str()), mt = std::atoi(line.substr(72, 3).c_str());
    if (mf && mt) out[{mf, mt}].push_back(line);
  }
  return out;
}

class Reader {
 public:
  explicit Reader(const Lines& lines) : lines_(lines) {}

  std::vector<double> cont() {
    const std::string& line = lines_.at(i_++);
    std::vector<double> v;
    for (int k = 0; k < 66; k += 11) v.push_back(number(line.substr(k, 11)));
    return v;
  }

  std::vector<double> values(int n) {
    std::vector<double> v;
    while (int(v.size()) < n) {
      auto row = cont();
      v.insert(v.end(), row.begin(), row.end());
    }
    v.resize(n);
    return v;
  }

  // Interpolation laws are skipped: the data are used lin-lin throughout.
  std::vector<double> tab1(Points& pts) {
    auto c = cont();
    values(2 * int(c[4]));
    auto xy = values(2 * int(c[5]));
    pts.clear();
    for (std::size_t k = 0; k + 1 < xy.size(); k += 2) pts.push_back({xy[k], xy[k + 1]});
    return c;
  }

  std::vector<double> tab2() {
    auto c = cont();
    values(2 * int(c[4]));
    return c;
  }

  std::vector<double> list(std::vector<double>& items) {
    auto c = cont();
    items = values(int(c[4]));
    return c;
  }

 private:
  const Lines& lines_;
  std::size_t i_ = 0;
};

// P_l(x) for l = 0..lmax at each x, by the three-term recursion.
std::vector<std::vector<double>> legendreTable(int lmax, const std::vector<double>& xs) {
  std::vector<std::vector<double>> t{std::vector<double>(xs.size(), 1.0), xs};
  for (int l = 2; l <= lmax; ++l) {
    std::vector<double> row(xs.size());
    for (std::size_t i = 0; i < xs.size(); ++i)
      row[i] = ((2 * l - 1) * xs[i] * t[l - 1][i] - (l - 1) * t[l - 2][i]) / l;
    t.push_back(std::move(row));
  }
  return t;
}

// Edges of kBins equal-probability bins of a tabulated pdf(mu).
std::vector<double> equiprobable(const std::vector<double>& mu, std::vector<double> pdf) {
  for (double& p : pdf) p = std::max(p, 0.0);
  std::vector<double> cdf{0.0};
  for (std::size_t k = 1; k < mu.size(); ++k)
    cdf.push_back(cdf.back() + 0.5 * (pdf[k] + pdf[k - 1]) * (mu[k] - mu[k - 1]));
  double total = cdf.back();
  std::vector<double> edges{-1.0};
  std::size_t k = 0;
  for (int b = 1; b < kBins; ++b) {
    double target = total * b / kBins;
    while (cdf[k + 1] < target) ++k;
    double frac = cdf[k + 1] > cdf[k] ? (target - cdf[k]) / (cdf[k + 1] - cdf[k]) : 0.0;
    edges.push_back(mu[k] + frac * (mu[k + 1] - mu[k]));
  }
  edges.push_back(1.0);
  return edges;
}

std::vector<std::pair<double, std::vector<double>>> angular(const Lines& lines) {
  Reader r(lines);
  int ltt = int(r.cont()[3]);
  r.cont();  // LI, LCT (2 = centre of mass)
  std::vector<double> grid(401);
  for (int k = 0; k <= 400; ++k) grid[k] = -1 + 2.0 * k / 400;
  auto poly = legendreTable(100, grid);
  std::vector<std::pair<double, std::vector<double>>> table;
  if (ltt == 1 || ltt == 3) {
    int ne = int(r.tab2()[5]);
    for (int e = 0; e < ne; ++e) {
      std::vector<double> coeffs;
      auto c = r.list(coeffs);
      std::vector<double> pdf(grid.size(), 0.5);
      for (std::size_t l = 1; l <= coeffs.size(); ++l) {
        double w = (2.0 * l + 1) / 2 * coeffs[l - 1];
        for (std::size_t i = 0; i < grid.size(); ++i) pdf[i] = pdf[i] + w * poly[l][i];
      }
      table.push_back({c[1], equiprobable(grid, pdf)});
    }
  }
  if (ltt == 2 || ltt == 3) {
    int ne = int(r.tab2()[5]);
    for (int e = 0; e < ne; ++e) {
      Points pts;
      auto c = r.tab1(pts);
      std::vector<double> mu, pdf;
      for (auto [x, y] : pts) mu.push_back(x), pdf.push_back(y);
      table.push_back({c[1], equiprobable(mu, pdf)});
    }
  }
  return table;
}

// Points inside [kMinEnergy, kMaxEnergy], keeping one neighbour on each side.
Points clip(const Points& pts) {
  int first = -1, last = -1;
  for (int k = 0; k < int(pts.size()); ++k)
    if (pts[k].first >= kMinEnergy && pts[k].first <= kMaxEnergy) {
      if (first < 0) first = k;
      last = k;
    }
  if (first < 0) return {};
  int lo = std::max(first - 1, 0), hi = std::min(last + 1, int(pts.size()) - 1);
  return Points(pts.begin() + lo, pts.begin() + hi + 1);
}

std::string format(const char* f, double a, double b = 0) {
  char buf[64];
  std::snprintf(buf, sizeof buf, f, a, b);
  return buf;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::fprintf(stderr, "usage: extract_endf <evaluation.endf> <output.xs>\n");
    return 2;
  }
  auto sec = sections(argv[1]);
  for (auto key : {std::pair{3, 1}, {1, 452}, {4, 2}})
    if (!sec.count(key)) {
      std::fprintf(stderr, "extract_endf: MF%d MT%d missing\n", key.first, key.second);
      return 1;
    }
  std::string src = argv[1];
  std::vector<std::string> out{
      "# ENDF/B-VIII.0 fast-neutron data extracted by tools/extract_endf from " +
          src.substr(src.find_last_of('/') + 1),
      format("AWR %.6f", Reader(sec[{3, 1}]).cont()[1])};

  std::vector<int> wanted{2, 16, 18, 102};
  for (int mt = 51; mt <= 91; ++mt) wanted.push_back(mt);  // inelastic levels + continuum
  for (int mt : wanted) {
    if (!sec.count({3, mt})) continue;
    Reader r(sec[{3, mt}]);
    r.cont();
    Points pts;
    auto c = r.tab1(pts);
    pts = clip(pts);
    double peak = 0;
    for (auto [e, y] : pts) peak = std::max(peak, y);
    if (pts.empty() || peak <= 0) continue;
    char head[96];
    std::snprintf(head, sizeof head, "XS %d %.6e %zu", mt, c[1], pts.size());
    out.push_back(head);
    for (auto [e, y] : pts) out.push_back(format("%.6e %.6e", e, y));
  }

  Reader nu(sec[{1, 452}]);
  Points pts;
  if (int(nu.cont()[3]) == 2) {
    nu.tab1(pts);
  } else {  // polynomial in E
    std::vector<double> coeffs;
    nu.list(coeffs);
    for (int k = 0; k <= 50; ++k) {
      double e = kMinEnergy * std::pow(kMaxEnergy / kMinEnergy, k / 50.0), v = 0;
      for (std::size_t j = 0; j < coeffs.size(); ++j) v += coeffs[j] * std::pow(e, double(j));
      pts.push_back({e, v});
    }
  }
  pts = clip(pts);
  out.push_back("NU " + std::to_string(pts.size()));
  for (auto [e, v] : pts) out.push_back(format("%.6e %.6f", e, v));

  auto rows = angular(sec[{4, 2}]);
  rows.erase(std::remove_if(rows.begin(), rows.end(),
                            [](const auto& row) {
                              return row.first < kMinEnergy / 2 || row.first > kMaxEnergy * 1.5;
                            }),
             rows.end());
  out.push_back("ANG " + std::to_string(rows.size()) + " " + std::to_string(kBins + 1));
  for (const auto& [e, edges] : rows) {
    std::string line = format("%.6e", e);
    for (double x : edges) line += format(" %.5f", x);
    out.push_back(line);
  }

  std::FILE* f = std::fopen(argv[2], "w");
  if (!f) {
    std::fprintf(stderr, "extract_endf: cannot write %s\n", argv[2]);
    return 1;
  }
  for (const auto& line : out) std::fprintf(f, "%s\n", line.c_str());
  std::fclose(f);
  return 0;
}
