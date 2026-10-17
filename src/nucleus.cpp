#include "nucleus.hpp"

#include "constants.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <unordered_map>

namespace {

// SEMF coefficients, MeV (Rohlf, "Modern Physics", 1994).
constexpr double kVolume = 15.75;
constexpr double kSurface = 17.80;
constexpr double kCoulombTerm = 0.711;
constexpr double kAsymmetry = 23.70;
constexpr double kPairing = 11.18;

// Myers-Swiatecki (1967) liquid-drop surface parameters.
constexpr double kMSSurface = 17.9439;
constexpr double kMSKappa = 1.7826;

// Myers & Swiatecki, Nucl. Phys. 81, 1 (1966).
constexpr double kMSVolume = 15.677;
constexpr double kMSSurface66 = 18.56;
constexpr double kMSSymmetry = 1.79;
constexpr double kMSCoulomb = 0.717;
constexpr double kMSDiffuseness = 1.21129;
constexpr double kMSPairing = 11.0;
constexpr double kMSShellC = 5.8;
constexpr double kMSShellSmall = 0.325;

// F(n): how far n nucleons sit from a closed shell. Zero at magic numbers,
// largest mid-shell (the gap between a linear interpolation of 3/5 n^(5/3)
// across the shell and the function itself).
double shellF(int n) {
  static const int kMagic[] = {0, 2, 8, 14, 28, 50, 82, 126, 184};
  auto p53 = [](double x) { return std::pow(x, 5.0 / 3.0); };
  for (size_t i = 1; i < sizeof(kMagic) / sizeof(kMagic[0]); ++i) {
    if (n > kMagic[i]) continue;
    double lo = kMagic[i - 1], hi = kMagic[i];
    double q = 0.6 * (p53(hi) - p53(lo)) / (hi - lo);
    return q * (n - lo) - 0.6 * (p53(n) - p53(lo));
  }
  return 0.0;
}

struct MassEntry {
  double excessMeV;
  bool extrapolated;  // marked '#' in AME2020: from systematics, not measured
};

// Reads a fixed-width numeric field; false if it holds no number.
bool field(const std::string& line, size_t pos, size_t len, double& out) {
  if (line.size() < pos + len) return false;
  std::string text = line.substr(pos, len);
  for (char& c : text)
    if (c == '#') c = '.';
  char* end = nullptr;
  out = std::strtod(text.c_str(), &end);
  if (end == text.c_str()) return false;
  for (; *end; ++end)
    if (*end != ' ') return false;
  return true;
}

// The full AME2020 table (Wang et al., Chinese Phys. C 45, 030003), keyed by
// 1000 Z + A. Columns follow the Fortran format given in the file header.
const std::unordered_map<int, MassEntry>& massTable() {
  static const std::unordered_map<int, MassEntry> table = [] {
    const char* path = ATOMSIM_DATA_DIR "/mass_1.mas20";
    std::ifstream in(path);
    if (!in) {
      std::fprintf(stderr, "atomsim: cannot open AME2020 mass table %s\n", path);
      std::exit(1);
    }
    std::unordered_map<int, MassEntry> t;
    std::string line;
    double Z, A, excess;
    while (std::getline(in, line))
      if (field(line, 9, 5, Z) && field(line, 14, 5, A) && field(line, 28, 14, excess))
        t[1000 * int(Z) + int(A)] = {excess / 1000.0, line.substr(28, 14).find('#') != std::string::npos};
    return t;
  }();
  return table;
}


}  // namespace

double semfBindingEnergy(int Z, int A) {
  if (A <= 0) return 0.0;
  int N = A - Z;
  double a3 = std::cbrt(static_cast<double>(A));
  double be = kVolume * A - kSurface * a3 * a3 - kCoulombTerm * Z * (Z - 1) / a3 -
              kAsymmetry * double(N - Z) * (N - Z) / A;
  if (A % 2 == 0) be += (Z % 2 == 0 ? 1 : -1) * kPairing / std::sqrt(double(A));
  return be;
}

double shellBindingEnergy(int Z, int A) {
  int N = A - Z;
  double I = double(N - Z) / A, s = 1 - kMSSymmetry * I * I, a3 = std::cbrt(double(A));
  double energy = -kMSVolume * s * A + kMSSurface66 * s * a3 * a3 + kMSCoulomb * Z * Z / a3 -
                  kMSDiffuseness * double(Z) * Z / A;
  if (A % 2 == 0) energy += (Z % 2 == 0 ? -1 : 1) * kMSPairing / std::sqrt(double(A));
  energy += kMSShellC * ((shellF(N) + shellF(Z)) / std::pow(A / 2.0, 2.0 / 3.0) -
                         kMSShellSmall * a3);
  return -energy;
}

std::optional<double> measuredMassExcess(int Z, int A) {
  auto it = massTable().find(1000 * Z + A);
  if (it == massTable().end()) return std::nullopt;
  return it->second.excessMeV;
}

bool isMeasured(int Z, int A) {
  auto it = massTable().find(1000 * Z + A);
  return it != massTable().end() && !it->second.extrapolated;
}

double bindingEnergy(int Z, int A) {
  if (auto excess = measuredMassExcess(Z, A))
    return Z * phys::kExcessHydrogen + (A - Z) * phys::kExcessNeutron - *excess;
  return shellBindingEnergy(Z, A);
}

double neutronSeparation(int Z, int A) { return bindingEnergy(Z, A) - bindingEnergy(Z, A - 1); }

double betaQ(int Z, int A) {
  return phys::kNeutronMinusHydrogen + bindingEnergy(Z + 1, A) - bindingEnergy(Z, A);
}

double nuclearMass(int Z, int A) {
  return Z * phys::kMassProton + (A - Z) * phys::kMassNeutron - bindingEnergy(Z, A);
}

double nuclearRadius(int A) { return 1.2 * std::cbrt(static_cast<double>(A)); }

double fissility(int Z, int A) {
  double I = double(A - 2 * Z) / A;
  double critical = 50.883 * (1 - kMSKappa * I * I);  // (Z^2/A) at which the barrier vanishes
  return double(Z) * Z / A / critical;
}

double fissionBarrier(int Z, int A) {
  double I = double(A - 2 * Z) / A;
  double a3 = std::cbrt(static_cast<double>(A));
  double surface = kMSSurface * (1 - kMSKappa * I * I) * a3 * a3;
  double x = fissility(Z, A);
  if (x >= 1) return 0.0;
  if (x > 2.0 / 3.0) return 0.83 * std::pow(1 - x, 3) * surface;
  return 0.38 * (0.75 - x) * surface;
}

double hillWheeler(double excitation, double barrier, double hbarOmega) {
  return 1.0 / (1.0 + std::exp(2 * phys::kPi * (barrier - excitation) / hbarOmega));
}

const char* elementSymbol(int Z) {
  static const char* kSymbols[] = {
      "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si",
      "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
      "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru",
      "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
      "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
      "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
      "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm"};
  constexpr int kCount = sizeof(kSymbols) / sizeof(kSymbols[0]);
  return (Z >= 0 && Z < kCount) ? kSymbols[Z] : "?";
}
