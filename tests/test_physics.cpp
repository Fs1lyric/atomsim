// Physics checks: conservation laws, textbook values, and agreement with
// measured U-235(n_th,f) data. Each bound is the accuracy the model claims.
#include "constants.hpp"
#include "electrons.hpp"
#include "fission.hpp"
#include "hydrogen.hpp"
#include "nucleus.hpp"

#include <cmath>
#include <cstdio>
#include <random>
#include <string>

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...)                                                \
  do {                                                                  \
    ++checks;                                                           \
    if (!(cond)) {                                                      \
      ++failures;                                                       \
      std::printf("FAIL %s:%d  %s\n      ", __FILE__, __LINE__, #cond); \
      std::printf(__VA_ARGS__);                                         \
      std::printf("\n");                                                \
    }                                                                   \
  } while (0)

#define CHECK_NEAR(a, b, tol)                                                 \
  CHECK(std::fabs((a) - (b)) <= (tol), "%s = %.10g, expected %.10g +/- %.3g", \
        #a, (double)(a), (double)(b), (double)(tol))

static void testHydrogen() {
  const double mu = phys::kHydrogenReducedMass;  // in electron masses
  for (auto [n, l] : {std::pair{1, 0}, {2, 0}, {2, 1}, {3, 2}, {4, 1}}) {
    RadialState s = solveHydrogenic(n, l, 1.0, mu);
    double exact = schrodingerEnergyEV(n, 1.0, mu);
    CHECK_NEAR(s.energyEV, exact, 1e-5);
    // <r> = [3n^2 - l(l+1)] / (2 mu Z) Bohr radii
    CHECK_NEAR(s.meanRadius, (3.0 * n * n - l * (l + 1)) / (2.0 * mu), 1e-3);
  }
  CHECK_NEAR(schrodingerEnergyEV(1, 1.0, mu), -13.5983, 1e-3);
  // 2p3/2 - 2p1/2 Dirac splitting: mu c^2 alpha^4 / 32 to leading order
  // (the measured 10.969 GHz adds ~0.025 GHz of QED, which is not modelled).
  const double muMeV = phys::kHydrogenReducedMassMeV;
  double split = diracEnergyEV(2, 1.5, 1.0, muMeV) - diracEnergyEV(2, 0.5, 1.0, muMeV);
  double leading = muMeV * 1e6 * std::pow(phys::kAlpha, 4) / 32;
  CHECK(std::fabs(split / leading - 1) < 1e-4, "split %.6g eV vs %.6g eV", split, leading);
  CHECK(std::fabs(split / phys::kPlanckEVs / 1e9 - 10.969) < 0.04, "within QED of measured");
}

static void testNucleus() {
  CHECK_NEAR(bindingEnergy(26, 56), 492.26, 0.01);
  CHECK_NEAR(bindingEnergy(82, 208), 1636.43, 0.01);
  for (auto [Z, A] : {std::pair{26, 56}, {50, 120}, {82, 208}, {92, 236}}) {
    double rel = std::fabs(semfBindingEnergy(Z, A) / bindingEnergy(Z, A) - 1);
    CHECK(rel < 0.01, "SEMF off by %.2f%% for Z=%d A=%d", 100 * rel, Z, A);
  }
  // Fragment mass model vs AME2020 (doubly magic 132Sn, 140Xe, 96Sr, 208Pb).
  CHECK_NEAR(shellBindingEnergy(50, 132), 1102.85, 2.0);
  CHECK_NEAR(shellBindingEnergy(54, 140), 1160.7, 2.0);
  CHECK_NEAR(shellBindingEnergy(38, 96), 818.0, 2.0);
  CHECK_NEAR(shellBindingEnergy(82, 208), 1636.43, 2.0);
  CHECK_NEAR(neutronSeparation(92, 236), 6.5455, 0.002);
  CHECK_NEAR(neutronSeparation(92, 239), 4.8064, 0.002);
  CHECK_NEAR(neutronSeparation(94, 240), 6.534, 0.002);
  // Beta-stable end of the A=140 chain is Ce-140.
  CHECK(betaQ(56, 140) > 0 && betaQ(57, 140) > 0 && betaQ(58, 140) < 0,
        "A=140 chain should stop at Ce");
  CHECK(std::string(elementSymbol(92)) == "U", "symbol");
}

static void testBarrier() {
  Compound u235 = formCompound(92, 235, 0.0253e-6);
  Compound u238 = formCompound(92, 238, 0.0253e-6);
  Compound u238fast = formCompound(92, 238, 2.0);
  CHECK(u235.excitation > u235.barrier, "U-235 + thermal n must be above barrier");
  CHECK(u238.excitation < u238.barrier, "U-238 + thermal n must be below barrier");
  CHECK(u238.transmission < 1e-4, "U-238 thermal transmission %.3g", u238.transmission);
  CHECK(u238fast.transmission > 0.5, "U-238 + 2 MeV n should fission");
  CHECK(u235.barrier > 4.5 && u235.barrier < 7.0, "B_f(U-236) = %.2f", u235.barrier);
  CHECK_NEAR(hillWheeler(6.0, 6.0, 0.6), 0.5, 1e-12);
}

static void testElectrons() {
  CHECK(electronConfiguration(26) == "[Ar] 3d6 4s2", "%s", electronConfiguration(26).c_str());
  CHECK(electronConfiguration(92) == "[Rn] 5f3 6d1 7s2", "%s", electronConfiguration(92).c_str());
  CHECK(electronConfiguration(1) == "1s1", "%s", electronConfiguration(1).c_str());
}

static void testSingleEvent() {
  std::mt19937_64 rng(7);
  Compound cn = formCompound(92, 235, 0.0253e-6);
  for (int i = 0; i < 2000; ++i) {
    FissionEvent ev = simulateFission(cn, rng);
    CHECK(ev.light.Z + ev.heavy.Z == 92, "charge");
    CHECK(ev.light.A + ev.heavy.A == 236, "mass number");
    CHECK(ev.light.A <= ev.heavy.A, "ordering");
    CHECK_NEAR(ev.tke + ev.light.excitation + ev.heavy.excitation, ev.q + cn.excitation, 1e-9);
    // Back-to-back: equal and opposite momenta.
    double pL = ev.light.beta * (ev.light.mass + ev.light.kineticEnergy);
    double pH = ev.heavy.beta * (ev.heavy.mass + ev.heavy.kineticEnergy);
    CHECK_NEAR(pL, pH, 1e-6);
    CHECK_NEAR(ev.light.kineticEnergy + ev.heavy.kineticEnergy, ev.tke, 1e-9);
    for (const Fragment* f : {&ev.light, &ev.heavy}) {
      double spent = f->gammaEnergy;
      for (size_t k = 0; k < f->neutronCmEnergies.size(); ++k)
        spent += f->neutronCmEnergies[k] + f->neutronSeparations[k];
      CHECK_NEAR(spent, f->excitation, 1e-9);
      CHECK(f->finalA == f->A - (int)f->neutronCmEnergies.size(), "neutron count");
      // Rotational energy is locked into the gamma cascade; only the rest
      // can evaporate neutrons, and evaporation runs until it can't.
      CHECK(f->rotationalEnergy >= 0 && f->rotationalEnergy <= f->gammaEnergy + 1e-12,
            "rotational energy must leave as gammas");
      CHECK(f->gammaEnergy - f->rotationalEnergy < neutronSeparation(f->Z, f->finalA) + 1e-12,
            "evaporation stopped early");
      if (!f->betaChain.empty())
        CHECK(betaQ(f->betaChain.back().Z + 1, f->finalA) <= 0, "chain must end stable");
    }
    if (i == 0) {
      auto traj = coulombTrajectory(ev, 1e4);
      double e0 = traj.front().kineticLight + traj.front().kineticHeavy + traj.front().potential;
      for (auto& p : traj) {
        double e = p.kineticLight + p.kineticHeavy + p.potential;
        CHECK(std::fabs(e / e0 - 1) < 1e-7, "energy drift %.3g at t=%g", e / e0 - 1, p.time);
      }
      auto& end = traj.back();
      CHECK((end.kineticLight + end.kineticHeavy) / ev.tke > 0.99, "fragments not accelerated");
      CHECK_NEAR(end.kineticLight, ev.light.kineticEnergy, 0.01 * ev.tke);
    }
  }
}

// Averages against measured U-235(n_th,f) data:
//   <TKE> 170.5 MeV, <A_H> 139.6, prompt nu 2.414, <E_n,lab> ~2.0 MeV.
static void testAgainstExperiment() {
  std::mt19937_64 rng(11);
  Compound cn = formCompound(92, 235, 0.0253e-6);
  RunStats s = runMonteCarlo(cn, 20000, rng);
  CHECK(std::fabs(s.meanTke - 170.5) < 5.0, "<TKE> = %.2f", s.meanTke);
  CHECK(std::fabs(s.meanHeavyA - 139.6) < 1.0, "<A_H> = %.2f", s.meanHeavyA);
  CHECK(std::fabs(s.meanNu - 2.414) < 0.4, "<nu> = %.3f", s.meanNu);
  CHECK(std::fabs(s.meanNeutronLabEnergy - 2.0) < 0.4, "<E_n> = %.3f", s.meanNeutronLabEnergy);
  CHECK(std::fabs(s.meanPromptGamma - 6.75) < 1.0, "<E_gamma> = %.2f", s.meanPromptGamma);
  CHECK(std::fabs(s.meanSpin - 7.0) < 0.5, "<J> = %.2f", s.meanSpin);
}

int main() {
  testHydrogen();
  testNucleus();
  testBarrier();
  testElectrons();
  testSingleEvent();
  testAgainstExperiment();
  std::printf("%d/%d checks passed\n", checks - failures, checks);
  return failures == 0 ? 0 : 1;
}
