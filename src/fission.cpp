#include "fission.hpp"

#include "constants.hpp"
#include "nucleus.hpp"

#include <algorithm>
#include <cmath>

namespace {

// Curvature of the fission barrier for actinides, MeV (Bjørnholm & Lynn 1980).
constexpr double kBarrierCurvature = 0.6;

// Fission modes of the random-neck-rupture picture (Brosa et al., Phys. Rep.
// 197, 1990). The heavy-fragment peaks are pinned by shells (Standard I near
// doubly magic 132Sn, Standard II near N = 88) and barely move between
// actinides, so the same positions serve U-236, U-239 and Pu-240. Weights,
// means and widths are fits to U-235(n_th,f) yields.
struct Mode {
  double weight, heavyMean, sigma;
};
constexpr Mode kStandard1{0.25, 134.0, 2.6};
constexpr Mode kStandard2{0.748, 141.0, 4.6};
constexpr double kSuperlongWeight = 0.002;  // symmetric mode, centred on A/2
constexpr double kSuperlongSigma = 9.0;

// Charge: unchanged charge distribution shifted by charge polarisation; the
// heavy fragment carries ~0.5 fewer protons (Wahl systematics).
constexpr double kPolarisation = 0.5;
constexpr double kChargeWidth = 0.56;

// Effective centre-to-centre distance at scission in units of
// (A_L^1/3 + A_H^1/3) fm. Touching spheres would be 1.2; real scission shapes
// are elongated with a neck, giving ~1.76, the value at which the Coulomb
// energy reproduces Viola's TKE systematics (Phys. Rev. C 31, 1985).
constexpr double kScissionRadius = 1.76;
constexpr double kTkeWidth = 8.0;  // MeV, shape fluctuations at scission

double cbrtA(int A) { return std::cbrt(static_cast<double>(A)); }

int sampleHeavyMass(int A, std::mt19937_64& rng) {
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  double pick = uniform(rng), mean, sigma;
  if (pick < kSuperlongWeight) {
    mean = A / 2.0;
    sigma = kSuperlongSigma;
  } else if (pick < kSuperlongWeight + kStandard1.weight) {
    mean = kStandard1.heavyMean;
    sigma = kStandard1.sigma;
  } else {
    mean = kStandard2.heavyMean;
    sigma = kStandard2.sigma;
  }
  int heavy = int(std::lround(std::normal_distribution<double>(mean, sigma)(rng)));
  if (2 * heavy < A) heavy = A - heavy;
  return std::clamp(heavy, (A + 1) / 2, A - 60);
}

// Emits neutrons while the excitation exceeds the separation energy. Neutron
// energies in the fragment frame follow Weisskopf's evaporation spectrum
// eps * exp(-eps / T), with the residual-nucleus temperature from the Fermi-gas
// level density a = A / 11 MeV^-1 used for fission fragments by Madland & Nix
// (Nucl. Sci. Eng. 81, 1982). What remains goes out as gamma rays.
void evaporate(Fragment& f, std::mt19937_64& rng) {
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  double e = f.excitation;
  int A = f.A;
  for (;;) {
    double sn = neutronSeparation(f.Z, A);
    if (e <= sn) break;
    double maxEps = e - sn;
    double temperature = std::sqrt(maxEps / ((A - 1) / 11.0));
    double eps;
    do {  // Gamma(2, T) by the product of two uniforms, truncated at maxEps
      eps = -temperature * std::log((1.0 - uniform(rng)) * (1.0 - uniform(rng)));
    } while (eps > maxEps);

    // Isotropic in the fragment frame, then boosted by the fragment velocity.
    double v = std::sqrt(2 * eps / phys::kMassNeutron);  // units of c
    double cosT = 2 * uniform(rng) - 1;
    double vx = f.beta + v * cosT, vy = v * std::sqrt(1 - cosT * cosT);
    double labBeta2 = vx * vx + vy * vy;
    double lab = phys::kMassNeutron * (1 / std::sqrt(1 - labBeta2) - 1);

    f.neutronSeparations.push_back(sn);
    f.neutronCmEnergies.push_back(eps);
    f.neutronLabEnergies.push_back(lab);
    e -= sn + eps;
    --A;
  }
  f.gammaEnergy = e;
  f.finalA = A;
}

void betaDecayChain(Fragment& f) {
  for (int Z = f.Z; Z < 110; ++Z) {
    double q = betaQ(Z, f.finalA);
    if (q <= 0) break;
    f.betaChain.push_back({Z, q});
  }
}

}  // namespace

Compound formCompound(int targetZ, int targetA, double neutronEnergy) {
  Compound c;
  c.targetZ = targetZ;
  c.targetA = targetA;
  c.Z = targetZ;
  c.A = targetA + 1;
  c.neutronEnergy = neutronEnergy;
  c.separationEnergy = neutronSeparation(c.Z, c.A);
  // Only the centre-of-mass share of the neutron's energy can excite the nucleus.
  c.excitation = c.separationEnergy + neutronEnergy * targetA / (targetA + 1.0);
  c.fissility = fissility(c.Z, c.A);
  c.barrier = fissionBarrier(c.Z, c.A);
  c.transmission = hillWheeler(c.excitation, c.barrier, kBarrierCurvature);
  return c;
}

FissionEvent simulateFission(const Compound& cn, std::mt19937_64& rng) {
  FissionEvent ev;
  Fragment& L = ev.light;
  Fragment& H = ev.heavy;

  H.A = sampleHeavyMass(cn.A, rng);
  L.A = cn.A - H.A;
  double zMean = cn.Z * double(H.A) / cn.A - kPolarisation;
  H.Z = int(std::lround(std::normal_distribution<double>(zMean, kChargeWidth)(rng)));
  L.Z = cn.Z - H.Z;
  L.mass = nuclearMass(L.Z, L.A);
  H.mass = nuclearMass(H.Z, H.A);

  // Energy available: ground-state Q-value plus the compound excitation.
  ev.q = bindingEnergy(L.Z, L.A) + bindingEnergy(H.Z, H.A) - bindingEnergy(cn.Z, cn.A);
  double available = ev.q + cn.excitation;

  // TKE is the Coulomb energy of the two charges at scission.
  double coulomb = phys::kCoulomb * L.Z * H.Z / (kScissionRadius * (cbrtA(L.A) + cbrtA(H.A)));
  ev.tke = std::normal_distribution<double>(coulomb, kTkeWidth)(rng);
  ev.tke = std::min(std::max(ev.tke, 0.5 * coulomb), available - 0.1);
  ev.scissionDistance = phys::kCoulomb * L.Z * H.Z / ev.tke;

  // The rest is internal excitation, shared in proportion to the level-density
  // parameter (equal temperatures) of each fragment.
  double txe = available - ev.tke;
  L.excitation = txe * L.A / cn.A;
  H.excitation = txe - L.excitation;

  // Back-to-back with equal momenta: exact relativistic two-body split of TKE.
  double M = L.mass + H.mass, W = M + ev.tke;
  double p = std::sqrt(ev.tke * (ev.tke + 2 * M) * (2 * H.mass + ev.tke) * (2 * L.mass + ev.tke)) /
             (2 * W);
  for (Fragment* f : {&L, &H}) {
    double energy = std::hypot(p, f->mass);
    f->kineticEnergy = energy - f->mass;
    f->beta = p / energy;
  }
  // Fix rounding so the two kinetic energies sum to TKE exactly.
  H.kineticEnergy = ev.tke - L.kineticEnergy;

  evaporate(L, rng);
  evaporate(H, rng);
  betaDecayChain(L);
  betaDecayChain(H);
  return ev;
}

std::vector<TrajectoryPoint> coulombTrajectory(const FissionEvent& ev, double endTime) {
  const double mL = ev.light.mass, mH = ev.heavy.mass;
  const double k = phys::kCoulomb * ev.light.Z * ev.heavy.Z;  // MeV fm
  const double c = phys::kSpeedOfLight;

  // State: separation r (fm) and momentum P = |p| c (MeV) of each fragment in the
  // centre-of-mass frame. dP/dt = c k / r^2, dr/dt = c (beta_L + beta_H).
  auto rate = [&](double r, double P, double& dr, double& dP) {
    dr = c * (P / std::hypot(P, mL) + P / std::hypot(P, mH));
    dP = c * k / (r * r);
  };
  auto sample = [&](double t, double r, double P) {
    return TrajectoryPoint{t, r, std::hypot(P, mL) - mL, std::hypot(P, mH) - mH, k / r};
  };

  std::vector<double> times = {0};
  for (double t = 1e-3; t <= endTime * 1.0001; t *= std::sqrt(10.0)) times.push_back(t);

  double r = ev.scissionDistance, P = 0, t = 0;
  // Speed scale used to pick steps: fragments' final relative speed.
  double vScale = c * (ev.light.beta + ev.heavy.beta);
  std::vector<TrajectoryPoint> out;
  for (double target : times) {
    while (t < target) {
      double dt = std::min(2e-3 * r / vScale, target - t);
      double r1, p1, r2, p2, r3, p3, r4, p4;
      rate(r, P, r1, p1);
      rate(r + 0.5 * dt * r1, P + 0.5 * dt * p1, r2, p2);
      rate(r + 0.5 * dt * r2, P + 0.5 * dt * p2, r3, p3);
      rate(r + dt * r3, P + dt * p3, r4, p4);
      r += dt / 6 * (r1 + 2 * r2 + 2 * r3 + r4);
      P += dt / 6 * (p1 + 2 * p2 + 2 * p3 + p4);
      t += dt;
    }
    out.push_back(sample(target, r, P));
  }
  return out;
}

RunStats runMonteCarlo(const Compound& cn, int events, std::mt19937_64& rng) {
  RunStats s;
  s.events = events;
  long neutronCount = 0;
  for (int i = 0; i < events; ++i) {
    FissionEvent ev = simulateFission(cn, rng);
    s.meanQ += ev.q;
    s.meanTke += ev.tke;
    s.meanTxe += ev.light.excitation + ev.heavy.excitation;
    s.meanLightKe += ev.light.kineticEnergy;
    s.meanHeavyKe += ev.heavy.kineticEnergy;
    s.meanLightA += ev.light.A;
    s.meanHeavyA += ev.heavy.A;
    int nu = ev.neutrons();
    s.meanNu += nu;
    s.multiplicity[std::min<size_t>(nu, s.multiplicity.size() - 1)] += 1;
    for (const Fragment* f : {&ev.light, &ev.heavy}) {
      s.yieldPre[f->A] += 1;
      s.yieldPost[f->finalA] += 1;
      for (size_t n = 0; n < f->neutronLabEnergies.size(); ++n) {
        double e = f->neutronLabEnergies[n];
        s.meanNeutronLabEnergy += e;
        s.meanNeutronKe += f->neutronCmEnergies[n];
        s.meanNeutronBinding += f->neutronSeparations[n];
        s.neutronSpectrum[std::min<size_t>(size_t(e / 0.25), s.neutronSpectrum.size() - 1)] += 1;
        ++neutronCount;
      }
      s.meanPromptGamma += f->gammaEnergy;
      for (const auto& b : f->betaChain) s.meanBetaChain += b.q;
    }
  }
  double perEvent = 1.0 / events;
  for (double* v : {&s.meanQ, &s.meanTke, &s.meanTxe, &s.meanLightKe, &s.meanHeavyKe,
                    &s.meanLightA, &s.meanHeavyA, &s.meanNu, &s.meanNeutronKe,
                    &s.meanNeutronBinding, &s.meanPromptGamma, &s.meanBetaChain})
    *v *= perEvent;
  for (double& y : s.yieldPre) y *= perEvent;
  for (double& y : s.yieldPost) y *= perEvent;
  for (double& m : s.multiplicity) m *= perEvent;
  if (neutronCount > 0) {
    s.meanNeutronLabEnergy /= neutronCount;
    for (double& n : s.neutronSpectrum) n /= neutronCount;
  }
  return s;
}
