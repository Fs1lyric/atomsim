#pragma once

#include <optional>

// Nuclear binding energies, in MeV, positive for bound nuclei.

// Bethe-Weizsäcker semi-empirical mass formula (liquid drop + pairing).
double semfBindingEnergy(int Z, int A);

// Myers-Swiatecki (1966) liquid drop with their analytic shell correction:
// the mass model used for fission fragments, which sit near the N = 82 and
// Z = 50 shells where the plain SEMF is off by up to 10 MeV.
double shellBindingEnergy(int Z, int A);

// Measured mass excess from AME2020, for the handful of nuclides tabulated here.
std::optional<double> measuredMassExcess(int Z, int A);

// Measured binding energy when tabulated, otherwise shellBindingEnergy.
double bindingEnergy(int Z, int A);
bool isMeasured(int Z, int A);

// Neutron separation energy S_n of nucleus (Z, A): energy to remove one neutron.
double neutronSeparation(int Z, int A);

// Atomic-mass Q-value for beta-minus decay (Z, A) -> (Z+1, A).
double betaQ(int Z, int A);

// Nuclear (not atomic) rest mass, MeV/c^2.
double nuclearMass(int Z, int A);

// Sharp-surface radius 1.2 A^(1/3) fm.
double nuclearRadius(int A);

// Liquid-drop fissility x = E_Coulomb / (2 E_surface) of the spherical drop,
// with the Myers-Swiatecki (1967) surface-symmetry correction.
double fissility(int Z, int A);

// Liquid-drop fission barrier from the Cohen-Swiatecki (1963) fit to exact
// saddle-point shapes, scaled by the Myers-Swiatecki surface energy.
double fissionBarrier(int Z, int A);

// Hill-Wheeler transmission through an inverted-parabola barrier of height
// barrier and curvature hbarOmega: 1 / (1 + exp(2 pi (B - E) / hbar omega)).
double hillWheeler(double excitation, double barrier, double hbarOmega);

const char* elementSymbol(int Z);
