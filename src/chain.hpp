#pragma once

#include <cstdint>
#include <random>
#include <vector>

// Neutron chain reaction in a bare sphere of fissile metal: Monte Carlo
// transport of fast neutrons with one-group (fission-spectrum averaged)
// cross-sections.

struct Material {
  double density;        // g/cm^3
  double molarMass;      // g/mol
  double sigmaFission;   // barns
  double sigmaScatter;   // elastic, barns; isotropic in the lab frame
  double sigmaCapture;   // barns
  double nu;             // mean neutrons per fission
  double neutronSpeed;   // cm/s
};

// Pure U-235 metal. Cross-sections from B. C. Reed, "The Physics of the
// Manhattan Project" (Springer, 4th ed.), for ~2 MeV fission neutrons.
Material uranium235();

double numberDensity(const Material& m);             // nuclei per cm^3
double sphereRadius(const Material& m, double massKg);  // cm

// Effective multiplication factor of a bare sphere by the generation method.
double multiplicationFactor(const Material& m, double massKg, int neutronsPerGeneration,
                            std::uint64_t seed);

// Mass at which k = 1, by bisection.
double criticalMass(const Material& m, std::uint64_t seed);

struct Neutron {
  double x, y, z;     // cm from the centre
  double dx, dy, dz;  // unit direction
  double toCollision; // cm left to the next collision
};

struct Flash {        // a fission site, for drawing
  double x, y, z;
};

// Time-dependent chain reaction. When the tracked population exceeds the cap,
// half the neutrons are dropped and each survivor stands for twice as many.
class ChainReaction {
 public:
  ChainReaction(const Material& m, double massKg, std::uint64_t seed, int cap = 4000);

  void inject(int count);           // stray neutrons appearing at the centre
  void advance(double seconds);     // move every neutron forward in time

  const std::vector<Neutron>& neutrons() const { return neutrons_; }
  const std::vector<Flash>& recentFissions() const { return flashes_; }  // last advance
  double radius() const { return radius_; }
  double time() const { return time_; }               // s
  double population() const { return neutrons_.size() * weight_; }
  double fissions() const { return fissions_; }       // total, weighted

 private:
  Neutron spawn(double x, double y, double z);
  double distanceToSurface(const Neutron& n) const;

  Material material_;
  double radius_, sigmaTotal_, weight_ = 1, time_ = 0, fissions_ = 0;
  int cap_;
  std::mt19937_64 rng_;
  std::vector<Neutron> neutrons_;
  std::vector<Flash> flashes_;
};
