#include "nucleus.hpp"

#include "constants.hpp"

#include <cmath>

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

struct MassEntry {
  int Z, A;
  double excessKeV;
};

// AME2020 mass excesses (Wang et al., Chinese Phys. C 45, 030003).
constexpr MassEntry kMassTable[] = {
    {1, 1, 7288.971},     {2, 4, 2424.916},     {6, 12, 0.0},
    {26, 56, -60607.1},   {82, 208, -21748.6},  {92, 235, 40920.4},
    {92, 236, 42446.3},   {92, 238, 47308.9},   {92, 239, 50573.9},
    {94, 239, 48589.9},   {94, 240, 50127.0},
};

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

std::optional<double> measuredMassExcess(int Z, int A) {
  for (const auto& e : kMassTable)
    if (e.Z == Z && e.A == A) return e.excessKeV / 1000.0;
  return std::nullopt;
}

bool isMeasured(int Z, int A) { return measuredMassExcess(Z, A).has_value(); }

double bindingEnergy(int Z, int A) {
  if (auto excess = measuredMassExcess(Z, A))
    return Z * phys::kExcessHydrogen + (A - Z) * phys::kExcessNeutron - *excess;
  return semfBindingEnergy(Z, A);
}

double neutronSeparation(int Z, int A) { return bindingEnergy(Z, A) - bindingEnergy(Z, A - 1); }

double semfNeutronSeparation(int Z, int A) {
  return semfBindingEnergy(Z, A) - semfBindingEnergy(Z, A - 1);
}

double betaQ(int Z, int A) {
  return phys::kNeutronMinusHydrogen + semfBindingEnergy(Z + 1, A) - semfBindingEnergy(Z, A);
}

double nuclearMass(int Z, int A) {
  return Z * phys::kMassProton + (A - Z) * phys::kMassNeutron - semfBindingEnergy(Z, A);
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
