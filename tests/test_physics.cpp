// Physics checks: conservation laws, textbook values, and agreement with
// measured U-235(n_th,f) data. Each bound is the accuracy the model claims.
#include "chain.hpp"
#include "constants.hpp"
#include "electrons.hpp"
#include "fission.hpp"
#include "hydrogen.hpp"
#include "nucleus.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <random>
#include <string>

namespace {

// printf-style failure message.
__attribute__((format(printf, 1, 2))) std::string describe(const char* fmt, ...) {
  char buf[256];
  va_list args;
  va_start(args, fmt);
  std::vsnprintf(buf, sizeof buf, fmt, args);
  va_end(args);
  return buf;
}

}  // namespace

TEST(Physics, Hydrogen) {
  const double mu = phys::kHydrogenReducedMass;  // in electron masses
  for (auto [n, l] : {std::pair{1, 0}, {2, 0}, {2, 1}, {3, 2}, {4, 1}}) {
    RadialState s = solveHydrogenic(n, l, 1.0, mu);
    double exact = schrodingerEnergyEV(n, 1.0, mu);
    EXPECT_NEAR(s.energyEV, exact, 1e-5);
    // <r> = [3n^2 - l(l+1)] / (2 mu Z) Bohr radii
    EXPECT_NEAR(s.meanRadius, (3.0 * n * n - l * (l + 1)) / (2.0 * mu), 1e-3);
  }
  EXPECT_NEAR(schrodingerEnergyEV(1, 1.0, mu), -13.5983, 1e-3);
  // 2p3/2 - 2p1/2 Dirac splitting: mu c^2 alpha^4 / 32 to leading order
  // (the measured 10.969 GHz adds ~0.025 GHz of QED, which is not modelled).
  const double muMeV = phys::kHydrogenReducedMassMeV;
  double split = diracEnergyEV(2, 1.5, 1.0, muMeV) - diracEnergyEV(2, 0.5, 1.0, muMeV);
  double leading = muMeV * 1e6 * std::pow(phys::kAlpha, 4) / 32;
  EXPECT_TRUE(std::fabs(split / leading - 1) < 1e-4) << describe("split %.6g eV vs %.6g eV", split, leading);
  EXPECT_TRUE(std::fabs(split / phys::kPlanckEVs / 1e9 - 10.969) < 0.04) << "within QED of measured";
}

TEST(Physics, Nucleus) {
  EXPECT_NEAR(bindingEnergy(26, 56), 492.26, 0.01);
  // Fragments come from the full AME2020 table, not the model.
  EXPECT_TRUE(isMeasured(54, 140) && isMeasured(38, 96) && isMeasured(50, 132)) << "AME2020 loaded";
  // AME prints B/A to 0.1 eV, so B for A = 140 is only known to 7 eV from that column.
  EXPECT_NEAR(bindingEnergy(54, 140), 8.2908875 * 140, 1e-5);
  EXPECT_NEAR(bindingEnergy(82, 208), 1636.43, 0.01);
  for (auto [Z, A] : {std::pair{26, 56}, {50, 120}, {82, 208}, {92, 236}}) {
    double rel = std::fabs(semfBindingEnergy(Z, A) / bindingEnergy(Z, A) - 1);
    EXPECT_TRUE(rel < 0.01) << describe("SEMF off by %.2f%% for Z=%d A=%d", 100 * rel, Z, A);
  }
  // Fragment mass model vs AME2020 (doubly magic 132Sn, 140Xe, 96Sr, 208Pb).
  EXPECT_NEAR(shellBindingEnergy(50, 132), 1102.85, 2.0);
  EXPECT_NEAR(shellBindingEnergy(54, 140), 1160.7, 2.0);
  EXPECT_NEAR(shellBindingEnergy(38, 96), 818.0, 2.0);
  EXPECT_NEAR(shellBindingEnergy(82, 208), 1636.43, 2.0);
  EXPECT_NEAR(neutronSeparation(92, 236), 6.5455, 0.002);
  EXPECT_NEAR(neutronSeparation(92, 239), 4.8064, 0.002);
  EXPECT_NEAR(neutronSeparation(94, 240), 6.5342, 0.002);
  // Beta-stable end of the A=140 chain is Ce-140.
  EXPECT_TRUE(betaQ(56, 140) > 0 && betaQ(57, 140) > 0 && betaQ(58, 140) < 0) << "A=140 chain should stop at Ce";
  EXPECT_TRUE(std::string(elementSymbol(92)) == "U") << "symbol";
}

TEST(Physics, Barrier) {
  Compound u235 = formCompound(92, 235, 0.0253e-6);
  Compound u238 = formCompound(92, 238, 0.0253e-6);
  Compound u238fast = formCompound(92, 238, 2.0);
  EXPECT_TRUE(u235.excitation > u235.barrier) << "U-235 + thermal n must be above barrier";
  EXPECT_TRUE(u238.excitation < u238.barrier) << "U-238 + thermal n must be below barrier";
  EXPECT_TRUE(u238.transmission < 1e-4) << describe("U-238 thermal transmission %.3g", u238.transmission);
  EXPECT_TRUE(u238fast.transmission > 0.5) << "U-238 + 2 MeV n should fission";
  EXPECT_TRUE(u235.barrier > 4.5 && u235.barrier < 7.0) << describe("B_f(U-236) = %.2f", u235.barrier);
  EXPECT_NEAR(hillWheeler(6.0, 6.0, 0.6), 0.5, 1e-12);
}

TEST(Physics, Electrons) {
  EXPECT_TRUE(electronConfiguration(26) == "[Ar] 3d6 4s2") << describe("%s", electronConfiguration(26).c_str());
  EXPECT_TRUE(electronConfiguration(92) == "[Rn] 5f3 6d1 7s2") << describe("%s", electronConfiguration(92).c_str());
  EXPECT_TRUE(electronConfiguration(1) == "1s1") << describe("%s", electronConfiguration(1).c_str());
}

TEST(Physics, SingleEvent) {
  std::mt19937_64 rng(7);
  Compound cn = formCompound(92, 235, 0.0253e-6);
  for (int i = 0; i < 2000; ++i) {
    FissionEvent ev = simulateFission(cn, rng);
    EXPECT_TRUE(ev.light.Z + ev.heavy.Z == 92) << "charge";
    EXPECT_TRUE(ev.light.A + ev.heavy.A == 236) << "mass number";
    EXPECT_TRUE(ev.light.A <= ev.heavy.A) << "ordering";
    EXPECT_NEAR(ev.tke + ev.light.excitation + ev.heavy.excitation, ev.q + cn.excitation, 1e-9);
    // Back-to-back: equal and opposite momenta.
    double pL = ev.light.beta * (ev.light.mass + ev.light.kineticEnergy);
    double pH = ev.heavy.beta * (ev.heavy.mass + ev.heavy.kineticEnergy);
    EXPECT_NEAR(pL, pH, 1e-6);
    EXPECT_NEAR(ev.light.kineticEnergy + ev.heavy.kineticEnergy, ev.tke, 1e-9);
    for (const Fragment* f : {&ev.light, &ev.heavy}) {
      double spent = f->gammaEnergy;
      for (size_t k = 0; k < f->neutronCmEnergies.size(); ++k)
        spent += f->neutronCmEnergies[k] + f->neutronSeparations[k];
      EXPECT_NEAR(spent, f->excitation, 1e-9);
      EXPECT_TRUE(f->finalA == f->A - (int)f->neutronCmEnergies.size()) << "neutron count";
      // Rotational energy is locked into the gamma cascade; only the rest
      // can evaporate neutrons, and evaporation runs until it can't.
      EXPECT_TRUE(f->rotationalEnergy >= 0 && f->rotationalEnergy <= f->gammaEnergy + 1e-12) << "rotational energy must leave as gammas";
      EXPECT_TRUE(f->gammaEnergy - f->rotationalEnergy < neutronSeparation(f->Z, f->finalA) + 1e-12) << "evaporation stopped early";
      if (!f->betaChain.empty()) {
        EXPECT_TRUE(betaQ(f->betaChain.back().Z + 1, f->finalA) <= 0) << "chain must end stable";
      }
    }
    if (i == 0) {
      auto traj = coulombTrajectory(ev, 1e4);
      double e0 = traj.front().kineticLight + traj.front().kineticHeavy + traj.front().potential;
      for (auto& p : traj) {
        double e = p.kineticLight + p.kineticHeavy + p.potential;
        EXPECT_TRUE(std::fabs(e / e0 - 1) < 1e-7) << describe("energy drift %.3g at t=%g", e / e0 - 1, p.time);
      }
      auto& end = traj.back();
      EXPECT_TRUE((end.kineticLight + end.kineticHeavy) / ev.tke > 0.99) << "fragments not accelerated";
      EXPECT_NEAR(end.kineticLight, ev.light.kineticEnergy, 0.01 * ev.tke);
    }
  }
}

// Averages against measured U-235(n_th,f) data:
//   <TKE> 170.5 MeV, <A_H> 139.6, prompt nu 2.414, <E_n,lab> ~2.0 MeV.
TEST(Physics, AgainstExperiment) {
  Compound cn = formCompound(92, 235, 0.0253e-6);
  RunStats s = runMonteCarlo(cn, 20000, 11);
  // Every averaged observable within 5% of experiment.
  auto within5 = [](double model, double measured) {
    return std::fabs(model / measured - 1) < 0.05;
  };
  EXPECT_TRUE(within5(s.meanTke, 170.5)) << describe("<TKE> = %.2f", s.meanTke);
  EXPECT_TRUE(within5(s.meanHeavyA, 139.6)) << describe("<A_H> = %.2f", s.meanHeavyA);
  EXPECT_TRUE(within5(s.meanLightA, 96.4)) << describe("<A_L> = %.2f", s.meanLightA);
  EXPECT_TRUE(within5(s.meanNu, 2.414)) << describe("<nu> = %.3f", s.meanNu);
  EXPECT_TRUE(within5(s.meanNeutronLabEnergy, 2.0)) << describe("<E_n> = %.3f", s.meanNeutronLabEnergy);
  EXPECT_TRUE(within5(s.meanPromptGamma, 6.75)) << describe("<E_gamma> = %.2f", s.meanPromptGamma);
  EXPECT_TRUE(within5(s.meanSpin, 7.0)) << describe("<J> = %.2f", s.meanSpin);
}

TEST(Physics, MonteCarloIsDeterministicAcrossThreadCounts) {
  Compound cn = formCompound(92, 235, 0.0253e-6);
  RunStats one = runMonteCarlo(cn, 5000, 3, 1);
  RunStats many = runMonteCarlo(cn, 5000, 3, 4);
  EXPECT_EQ(one.meanNu, many.meanNu);
  EXPECT_EQ(one.meanTke, many.meanTke);
  EXPECT_EQ(one.neutronCount, many.neutronCount);
  EXPECT_EQ(one.yieldPre, many.yieldPre);
  EXPECT_NE(runMonteCarlo(cn, 5000, 4, 4).meanTke, one.meanTke) << "seed must matter";
}

// Godiva (ICSBEP HEU-MET-FAST-001): a bare 52.42 kg sphere of 93.7% U-235
// metal at 18.74 g/cm^3 is exactly critical.
TEST(ChainReaction, ReproducesGodivaBenchmark) {
  EXPECT_NEAR(neutronSpeed(2.0), 1.95e9, 0.01e9);
  Fuel g = godiva();
  EXPECT_NEAR(sphereRadius(g, 52.42), 8.7407, 0.01);
  double k = multiplicationFactor(g, 52.42, 20000, 1);
  double mass = criticalMass(g, 2);
  std::printf("[ info ] Godiva k_eff %.4f (benchmark 1.0000), critical mass %.2f kg (52.42)\n",
              k, mass);
  EXPECT_NEAR(k, 1.0, 0.015);
  EXPECT_NEAR(mass, 52.42, 0.05 * 52.42);
}

TEST(ChainReaction, NaturalUraniumNeverGoesCritical) {
  // Fast fission in U-238 needs > 1 MeV, and inelastic scattering drains that.
  EXPECT_LT(multiplicationFactor(uranium(0.0072), 2000, 2000, 3), 0.5);
}

TEST(ChainReaction, GrowsAboveCriticalAndDiesBelow) {
  Fuel g = godiva();
  ChainReaction hot(g, 90, 4), cold(g, 20, 5);
  hot.inject(500);
  cold.inject(500);
  for (int i = 0; i < 300; ++i) {  // 300 ns
    hot.advance(1e-9);
    cold.advance(1e-9);
  }
  std::printf("[ info ] after 300 ns: 90 kg -> %.3g neutrons, 20 kg -> %.3g\n",
              hot.population(), cold.population());
  EXPECT_GT(hot.population(), 5000) << "90 kg must run away";
  EXPECT_LT(cold.population(), 50) << "20 kg must fizzle";
  EXPECT_GT(hot.fissions(), cold.fissions());
  for (const Neutron& n : hot.neutrons())
    EXPECT_LE(n.x * n.x + n.y * n.y + n.z * n.z, hot.radius() * hot.radius() * (1 + 1e-9));
}
