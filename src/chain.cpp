#include "chain.hpp"

#include "constants.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace {

constexpr double kAvogadro = 6.02214076e23;
constexpr double kSpeedOfLightCm = 2.99792458e10;

// Below 25 keV the extracted tables stop (resolved resonances need MF2
// reconstruction). Fast systems put ~1% of collisions there; they use the
// 25 keV cross-sections.
constexpr double kMinEnergy = 0.025;  // MeV

using Rng = std::mt19937_64;

double uniform(Rng& rng) { return std::uniform_real_distribution<double>(0.0, 1.0)(rng); }

void isotropic(Rng& rng, double& dx, double& dy, double& dz) {
  double cosT = 2 * uniform(rng) - 1, phi = 2 * phys::kPi * uniform(rng);
  double sinT = std::sqrt(1 - cosT * cosT);
  dx = sinT * std::cos(phi);
  dy = sinT * std::sin(phi);
  dz = cosT;
}

// Turns direction (u, v, w) through polar cosine mu and a random azimuth.
void rotate(Rng& rng, double& u, double& v, double& w, double mu) {
  double phi = 2 * phys::kPi * uniform(rng), s = std::sqrt(std::max(0.0, 1 - mu * mu));
  double c = std::cos(phi), sn = std::sin(phi);
  if (std::fabs(w) > 0.99999) {
    u = s * c;
    v = s * sn;
    w = w > 0 ? mu : -mu;
    return;
  }
  double t = std::sqrt(1 - w * w);
  double u2 = mu * u + s * (u * w * c - v * sn) / t;
  double v2 = mu * v + s * (v * w * c + u * sn) / t;
  w = mu * w - s * t * c;
  u = u2;
  v = v2;
}

int sampleCount(double mean, Rng& rng) {
  int base = int(mean);
  return base + (uniform(rng) < mean - base ? 1 : 0);
}

// Watt fission spectrum: a Maxwellian of temperature a, boosted (Romano et al.).
double sampleWatt(double a, double b, Rng& rng) {
  double c = std::cos(phys::kPi / 2 * uniform(rng));
  double w = -a * (std::log(1 - uniform(rng)) + std::log(1 - uniform(rng)) * c * c);
  return w + a * a * b / 4 + (2 * uniform(rng) - 1) * std::sqrt(a * a * b * w);
}

// Evaporation spectrum E exp(-E/theta), truncated at the available energy.
double sampleEvaporation(double theta, double maxEnergy, Rng& rng) {
  for (;;) {
    double e = -theta * std::log((1 - uniform(rng)) * (1 - uniform(rng)));
    if (e <= maxEnergy) return e;
  }
}

double exitDistance(const Neutron& n, double R) {
  double b = n.x * n.dx + n.y * n.dy + n.z * n.dz;
  double c = n.x * n.x + n.y * n.y + n.z * n.z - R * R;
  return -b + std::sqrt(std::max(0.0, b * b - c));
}

// Cross-sections are pre-tabulated on one logarithmic energy grid so a lookup
// is an index computation instead of a binary search per reaction.
constexpr double kGridMax = 20.0;  // MeV
constexpr int kGridPoints = 4000;

void gridPosition(double e, std::size_t& i, double& f) {
  double t = std::log(std::clamp(e, kMinEnergy, kGridMax) / kMinEnergy) /
             std::log(kGridMax / kMinEnergy) * (kGridPoints - 1);
  i = std::min<std::size_t>(std::size_t(t), kGridPoints - 2);
  f = t - double(i);
}

double gridEnergy(int i) {
  return kMinEnergy * std::pow(kGridMax / kMinEnergy, double(i) / (kGridPoints - 1));
}

double onGrid(const std::vector<double>& g, double e) {
  std::size_t i;
  double f;
  gridPosition(e, i, f);
  return g[i] + f * (g[i + 1] - g[i]);
}

double reactionSigma(const NuclideData& d, std::size_t r, double e) {
  return onGrid(d.reactionGrid[r], e);
}

double totalSigma(const NuclideData& d, double e) { return onGrid(d.totalGrid, e); }

double macroscopicTotal(const Fuel& m, double e) {  // 1/cm
  double s = 0;
  for (const auto& [data, density] : m.nuclides) s += density * totalSigma(*data, e);
  return s;
}

double sampleFlight(const Fuel& m, double e, Rng& rng) {
  return -std::log(1 - uniform(rng)) / macroscopicTotal(m, std::max(e, kMinEnergy));
}

Neutron born(double x, double y, double z, double energy, Rng& rng) {
  Neutron n{x, y, z, 0, 0, 1, energy, 0};
  isotropic(rng, n.dx, n.dy, n.dz);
  return n;
}

enum class Outcome { Moving, Fission, Absorbed };

// One collision. New neutrons (fission, (n,2n)) are appended to `out`.
Outcome collide(const Fuel& m, Neutron& n, Rng& rng, std::vector<Neutron>& out) {
  const double e = std::max(n.energy, kMinEnergy);
  double pick = uniform(rng) * macroscopicTotal(m, e);
  for (const auto& [dataPtr, density] : m.nuclides) {
    const NuclideData& d = *dataPtr;
    for (std::size_t r = 0; r < d.reactions.size(); ++r) {
      pick -= density * reactionSigma(d, r, e);
      if (pick > 0) continue;
      const Reaction& rx = d.reactions[r];
      const double A = d.awr;
      if (rx.mt == 18) {
        for (int k = sampleCount(d.nu.at(e), rng); k > 0; --k)
          out.push_back(born(n.x, n.y, n.z, sampleWatt(d.wattA, d.wattB, rng), rng));
        return Outcome::Fission;
      }
      if (rx.mt == 102) return Outcome::Absorbed;
      if (rx.mt == 2 || (rx.mt >= 51 && rx.mt <= 90)) {
        // Two-body kinematics on a target of A neutron masses, centre-of-mass
        // cosine mu; beta^2 = 1 + (A+1) Q / (A E) leaves the nucleus at level -Q.
        double mu = rx.mt == 2 ? d.sampleMuCm(e, rng) : 2 * uniform(rng) - 1;
        double beta2 = rx.mt == 2 ? 1.0 : 1 + (A + 1) * rx.q / (A * n.energy);
        double beta = std::sqrt(std::max(beta2, 0.0));
        double dd = 1 + A * A * beta * beta + 2 * A * beta * mu;
        double muLab = dd > 0 ? (1 + A * beta * mu) / std::sqrt(dd) : 0.0;
        n.energy *= dd / ((A + 1) * (A + 1));
        rotate(rng, n.dx, n.dy, n.dz, std::clamp(muLab, -1.0, 1.0));
        return Outcome::Moving;
      }
      // Continuum inelastic (91) and (n,2n) (16): evaporation from the excited
      // residual nucleus, Fermi-gas temperature with a = A/8 MeV^-1.
      double available = std::max(n.energy + rx.q * (A + 1) / A, 1e-3);
      double theta = std::sqrt(available / (A / 8));
      if (rx.mt == 16) {
        double first, second;
        do {
          first = sampleEvaporation(theta, available, rng);
          second = sampleEvaporation(theta, available, rng);
        } while (first + second > available);
        out.push_back(born(n.x, n.y, n.z, second, rng));
        n.energy = first;
      } else {
        n.energy = sampleEvaporation(theta, available, rng);
      }
      isotropic(rng, n.dx, n.dy, n.dz);
      return Outcome::Moving;
    }
  }
  return Outcome::Moving;  // rounding left pick > 0: treat as no interaction
}

[[noreturn]] void dataError(const std::string& path, const char* what) {
  std::fprintf(stderr, "atomsim: %s: %s\n", path.c_str(), what);
  std::exit(1);
}

}  // namespace

double Table::at(double e) const {
  if (x.empty() || e < x.front() || e > x.back()) return 0.0;
  auto it = std::upper_bound(x.begin(), x.end(), e);
  if (it == x.end()) return y.back();
  std::size_t i = std::size_t(it - x.begin());
  if (i == 0) return y.front();
  double f = (e - x[i - 1]) / (x[i] - x[i - 1]);
  return y[i - 1] + f * (y[i] - y[i - 1]);
}

NuclideData NuclideData::load(const std::string& path, double wattA, double wattB) {
  std::ifstream in(path);
  if (!in) dataError(path, "cannot open nuclear data");
  NuclideData d;
  d.wattA = wattA;
  d.wattB = wattB;
  std::string line, tag;
  auto readTable = [&](int n, Table& t) {
    for (int i = 0; i < n && std::getline(in, line); ++i) {
      double ev, v;
      std::istringstream(line) >> ev >> v;
      t.x.push_back(ev * 1e-6);  // eV -> MeV
      t.y.push_back(v);
    }
  };
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream head(line);
    head >> tag;
    if (tag == "AWR") {
      head >> d.awr;
    } else if (tag == "XS") {
      Reaction r;
      double qEv;
      int n;
      head >> r.mt >> qEv >> n;
      r.q = qEv * 1e-6;
      readTable(n, r.sigma);
      d.reactions.push_back(std::move(r));
    } else if (tag == "NU") {
      int n;
      head >> n;
      readTable(n, d.nu);
    } else if (tag == "ANG") {
      int n, edges;
      head >> n >> edges;
      if (edges != 33) dataError(path, "expected 33 angular bin edges");
      for (int i = 0; i < n && std::getline(in, line); ++i) {
        std::istringstream row(line);
        double ev;
        std::array<double, 33> bins{};
        row >> ev;
        for (double& b : bins) row >> b;
        d.angleEnergy.push_back(ev * 1e-6);
        d.angleBins.push_back(bins);
      }
    }
  }
  if (d.awr <= 0 || d.reactions.empty() || d.nu.x.empty() || d.angleBins.empty())
    dataError(path, "incomplete nuclear data");
  d.totalGrid.assign(kGridPoints, 0.0);
  for (const Reaction& r : d.reactions) {
    std::vector<double> g(kGridPoints);
    for (int i = 0; i < kGridPoints; ++i) {
      g[i] = r.sigma.at(gridEnergy(i));
      d.totalGrid[i] += g[i];
    }
    d.reactionGrid.push_back(std::move(g));
  }
  return d;
}

double NuclideData::sampleMuCm(double e, std::mt19937_64& rng) const {
  // Stochastic interpolation between the two tabulated incident energies.
  auto it = std::upper_bound(angleEnergy.begin(), angleEnergy.end(), e);
  std::size_t i = std::size_t(std::max<std::ptrdiff_t>(it - angleEnergy.begin(), 1)) - 1;
  if (i + 1 < angleEnergy.size()) {
    double f = (e - angleEnergy[i]) / (angleEnergy[i + 1] - angleEnergy[i]);
    if (uniform(rng) < f) ++i;
  }
  const auto& bins = angleBins[i];
  std::size_t b = std::min<std::size_t>(std::size_t(uniform(rng) * 32), 31);
  return bins[b] + uniform(rng) * (bins[b + 1] - bins[b]);
}

// Watt parameters: ENDF/B-VII fits at ~1 MeV incident energy.
const NuclideData& u235Data() {
  static const NuclideData d = NuclideData::load(ATOMSIM_DATA_DIR "/u235.xs", 0.988, 2.249);
  return d;
}
const NuclideData& u238Data() {
  static const NuclideData d = NuclideData::load(ATOMSIM_DATA_DIR "/u238.xs", 0.88111, 3.4005);
  return d;
}

Fuel uranium(double enrichment, double density) {
  double perGram = kAvogadro * 1e-24;  // atoms per barn-cm per (g/cm^3 / g/mol)
  return {density,
          {{&u235Data(), density * enrichment / 235.0439 * perGram},
           {&u238Data(), density * (1 - enrichment) / 238.0508 * perGram}}};
}

Fuel godiva() {
  // Benchmark atom densities (U-234, 4.9e-4 /b-cm, is left out).
  return {18.74, {{&u235Data(), 4.4994e-2}, {&u238Data(), 2.4984e-3}}};
}

double sphereRadius(const Fuel& m, double massKg) {
  return std::cbrt(3 * massKg * 1000 / (4 * phys::kPi * m.density));
}

double neutronSpeed(double energyMeV) {
  double gamma = 1 + energyMeV / phys::kMassNeutron;
  return kSpeedOfLightCm * std::sqrt(1 - 1 / (gamma * gamma));
}

double multiplicationFactor(const Fuel& m, double massKg, int perGeneration,
                            std::uint64_t seed) {
  const double R = sphereRadius(m, massKg);
  Rng rng(seed);
  const NuclideData& fuel = *m.nuclides.front().first;

  // Start uniformly in the sphere; the fission source converges in a few generations.
  std::vector<Neutron> source;
  for (int i = 0; i < perGeneration; ++i) {
    double r = R * std::cbrt(uniform(rng)), dx, dy, dz;
    isotropic(rng, dx, dy, dz);
    source.push_back(born(r * dx, r * dy, r * dz, sampleWatt(fuel.wattA, fuel.wattB, rng), rng));
  }

  constexpr int kInactive = 10, kActive = 40;
  double kSum = 0;
  std::vector<Neutron> bank, secondaries;
  for (int gen = 0; gen < kInactive + kActive; ++gen) {
    bank.clear();
    for (Neutron n : source) {
      std::vector<Neutron> stack{n};
      while (!stack.empty()) {
        Neutron cur = stack.back();
        stack.pop_back();
        for (;;) {
          double d = sampleFlight(m, cur.energy, rng);
          if (d >= exitDistance(cur, R)) break;  // leaks out
          cur.x += d * cur.dx;
          cur.y += d * cur.dy;
          cur.z += d * cur.dz;
          secondaries.clear();
          Outcome o = collide(m, cur, rng, secondaries);
          if (o == Outcome::Fission) {
            bank.insert(bank.end(), secondaries.begin(), secondaries.end());
            break;
          }
          stack.insert(stack.end(), secondaries.begin(), secondaries.end());  // (n,2n)
          if (o == Outcome::Absorbed) break;
        }
      }
    }
    double k = double(bank.size()) / source.size();
    if (gen >= kInactive) kSum += k;
    if (bank.empty()) return 0.0;
    for (auto& s : source) s = bank[std::size_t(uniform(rng) * bank.size()) % bank.size()];
  }
  return kSum / kActive;
}

double criticalMass(const Fuel& m, std::uint64_t seed) {
  double lo = 10, hi = 2000;
  for (int i = 0; i < 16; ++i) {
    double mid = std::sqrt(lo * hi);
    (multiplicationFactor(m, mid, 10000, seed + i) < 1 ? lo : hi) = mid;
  }
  return std::sqrt(lo * hi);
}

ChainReaction::ChainReaction(const Fuel& m, double massKg, std::uint64_t seed, int cap)
    : material_(m), radius_(sphereRadius(m, massKg)), cap_(cap), rng_(seed) {}

void ChainReaction::inject(int count) {
  const NuclideData& fuel = *material_.nuclides.front().first;
  for (int i = 0; i < count; ++i) {
    Neutron n = born(0, 0, 0, sampleWatt(fuel.wattA, fuel.wattB, rng_), rng_);
    n.toCollision = sampleFlight(material_, n.energy, rng_);
    neutrons_.push_back(n);
  }
}

void ChainReaction::advance(double seconds) {
  flashes_.clear();
  std::vector<Neutron> next, secondaries;
  next.reserve(neutrons_.size() * 2);
  for (Neutron n : neutrons_) {
    double timeLeft = seconds;
    bool alive = true;
    while (alive) {
      double reach = neutronSpeed(n.energy) * timeLeft;
      double step = std::min(reach, n.toCollision);
      if (exitDistance(n, radius_) <= step) {  // escapes the sphere
        alive = false;
        break;
      }
      n.x += step * n.dx;
      n.y += step * n.dy;
      n.z += step * n.dz;
      n.toCollision -= step;
      timeLeft -= step / neutronSpeed(n.energy);
      if (n.toCollision > 0) break;  // still in flight at the end of the step
      secondaries.clear();
      Outcome o = collide(material_, n, rng_, secondaries);
      for (Neutron s : secondaries) {
        s.toCollision = sampleFlight(material_, s.energy, rng_);
        next.push_back(s);
      }
      if (o == Outcome::Fission) {
        flashes_.push_back({n.x, n.y, n.z});
        fissions_ += weight_;
      }
      if (o != Outcome::Moving) alive = false;
      else n.toCollision = sampleFlight(material_, n.energy, rng_);
    }
    if (alive) next.push_back(n);
  }
  neutrons_.swap(next);
  while (int(neutrons_.size()) > cap_) {  // keep the work bounded; weights keep counts right
    std::shuffle(neutrons_.begin(), neutrons_.end(), rng_);
    neutrons_.resize(neutrons_.size() / 2);
    weight_ *= 2;
  }
  time_ += seconds;
}
