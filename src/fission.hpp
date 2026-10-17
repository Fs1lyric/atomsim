#pragma once

#include <array>
#include <cstdint>
#include <random>
#include <vector>

// Neutron-induced fission, from capture to beta-stable products.

struct Compound {
  int targetZ = 0, targetA = 0;
  int Z = 0, A = 0;               // compound nucleus (target + n)
  double neutronEnergy = 0;       // lab kinetic energy, MeV
  double separationEnergy = 0;    // S_n of the compound nucleus, MeV
  double excitation = 0;          // E* = S_n + centre-of-mass neutron energy, MeV
  double fissility = 0;
  double barrier = 0;             // liquid-drop B_f, MeV
  double transmission = 0;        // Hill-Wheeler barrier transmission
};

Compound formCompound(int targetZ, int targetA, double neutronEnergy);

struct BetaStep {
  int Z = 0;          // parent of this decay
  double q = 0;       // MeV, shared by electron, antineutrino and delayed gammas
};

struct Fragment {
  int Z = 0, A = 0;                       // at scission, before neutron emission
  double mass = 0;                        // nuclear rest mass, MeV/c^2
  double kineticEnergy = 0;               // after full Coulomb acceleration, MeV
  double beta = 0;                        // final speed / c
  double excitation = 0;                  // at scission, MeV
  int spin = 0;                           // angular momentum at scission, units of hbar
  double rotationalEnergy = 0;            // part of the excitation held by rotation, MeV
  std::vector<double> neutronSeparations; // S_n paid for each evaporated neutron
  std::vector<double> neutronCmEnergies;  // neutron energy in the fragment frame
  std::vector<double> neutronLabEnergies;
  double gammaEnergy = 0;                 // prompt gamma cascade, MeV
  int finalA = 0;                         // after prompt neutrons
  std::vector<BetaStep> betaChain;        // to the beta-stable isobar
};

struct FissionEvent {
  Fragment light, heavy;
  double q = 0;                // ground-state Q-value of CN -> light + heavy
  double tke = 0;              // total kinetic energy of the fragments
  double scissionDistance = 0; // centre-to-centre distance at scission, fm
  int neutrons() const {
    return int(light.neutronCmEnergies.size() + heavy.neutronCmEnergies.size());
  }
};

FissionEvent simulateFission(const Compound& cn, std::mt19937_64& rng);

// Relativistic two-body Coulomb acceleration after scission.
struct TrajectoryPoint {
  double time = 0;        // zs
  double separation = 0;  // fm
  double kineticLight = 0, kineticHeavy = 0, potential = 0;  // MeV
};
std::vector<TrajectoryPoint> coulombTrajectory(const FissionEvent& ev, double endTime);

struct RunStats {
  int events = 0;
  std::array<double, 260> yieldPre{};   // fragments per fission by A, before neutrons
  std::array<double, 260> yieldPost{};  // after prompt neutrons
  std::array<double, 10> multiplicity{};
  std::array<double, 40> neutronSpectrum{};  // lab energy, 0.25 MeV bins
  double meanQ = 0, meanTke = 0, meanTxe = 0, meanLightKe = 0, meanHeavyKe = 0;
  double meanLightA = 0, meanHeavyA = 0, meanNu = 0, meanNeutronLabEnergy = 0;
  double meanNeutronKe = 0, meanNeutronBinding = 0, meanPromptGamma = 0, meanBetaChain = 0;
  double meanSpin = 0;  // per fragment, hbar
  long neutronCount = 0;
};

// Runs `events` fissions on `threads` threads (0 = all cores). Deterministic in
// `seed`: the same seed gives bit-identical results for any thread count.
RunStats runMonteCarlo(const Compound& cn, int events, std::uint64_t seed, int threads = 0);
