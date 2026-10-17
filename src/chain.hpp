#pragma once

#include <array>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

// Neutron chain reaction in a bare sphere of uranium metal: continuous-energy
// Monte Carlo transport with ENDF/B-VIII.0 data (see tools/extract_endf.py).

// Linear-linear interpolated table; zero outside its energy range.
struct Table {
  std::vector<double> x, y;
  double at(double e) const;
};

struct Reaction {
  int mt;        // ENDF reaction number: 2 elastic, 16 (n,2n), 18 fission,
                 // 51-90 inelastic levels, 91 inelastic continuum, 102 capture
  double q;      // Q-value, MeV
  Table sigma;   // barns vs MeV
};

struct NuclideData {
  double awr = 0;                          // mass in neutron masses
  std::vector<Reaction> reactions;
  Table nu;                                // mean neutrons per fission
  std::vector<double> angleEnergy;         // MeV
  std::vector<std::array<double, 33>> angleBins;  // equiprobable mu_cm edges
  double wattA = 0, wattB = 0;             // fission spectrum, MeV and 1/MeV
  std::vector<double> totalGrid;           // sigma on the shared energy grid
  std::vector<std::vector<double>> reactionGrid;

  static NuclideData load(const std::string& path, double wattA, double wattB);
  double sampleMuCm(double e, std::mt19937_64& rng) const;
};

const NuclideData& u235Data();
const NuclideData& u238Data();

struct Fuel {
  double density;                          // g/cm^3
  std::vector<std::pair<const NuclideData*, double>> nuclides;  // atoms per barn-cm
};

// Uranium metal of the given U-235 weight fraction (U-238 for the rest).
Fuel uranium(double enrichment, double density = 18.74);
// Godiva, ICSBEP HEU-MET-FAST-001: 93.7 wt% U-235, critical at 52.42 kg.
Fuel godiva();

double sphereRadius(const Fuel& m, double massKg);  // cm
double neutronSpeed(double energyMeV);                  // cm/s

// k_eff of a bare sphere by the generation (power-iteration) method.
double multiplicationFactor(const Fuel& m, double massKg, int neutronsPerGeneration,
                            std::uint64_t seed);
// Mass at which k_eff = 1, by bisection.
double criticalMass(const Fuel& m, std::uint64_t seed);

struct Neutron {
  double x, y, z;      // cm from the centre
  double dx, dy, dz;   // unit direction
  double energy;       // MeV
  double toCollision;  // cm left to the next collision
};

struct Flash {         // a fission site, for drawing
  double x, y, z;
};

// Time-dependent chain reaction. When the tracked population exceeds the cap,
// half the neutrons are dropped and each survivor stands for twice as many.
class ChainReaction {
 public:
  ChainReaction(const Fuel& m, double massKg, std::uint64_t seed, int cap = 4000);

  void inject(int count);         // fission neutrons appearing at the centre
  void advance(double seconds);   // move every neutron forward in time

  const std::vector<Neutron>& neutrons() const { return neutrons_; }
  const std::vector<Flash>& recentFissions() const { return flashes_; }  // last advance
  double radius() const { return radius_; }
  double time() const { return time_; }                // s
  double population() const { return neutrons_.size() * weight_; }
  double fissions() const { return fissions_; }        // total, weighted

 private:
  Fuel material_;
  double radius_, weight_ = 1, time_ = 0, fissions_ = 0;
  int cap_;
  std::mt19937_64 rng_;
  std::vector<Neutron> neutrons_;
  std::vector<Flash> flashes_;
};
