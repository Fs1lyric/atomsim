#pragma once

#include <vector>

// A bound state of a one-electron atom, found by integrating the radial
// Schrödinger equation  u'' = [l(l+1)/r^2 - 2 mu Z / r - 2 mu E] u
// (atomic units, u = r R(r)) with Numerov's method and shooting on E.
struct RadialState {
  int n = 0, l = 0;
  double energyEV = 0;           // binding energy, eV (negative)
  double meanRadius = 0;         // <r> in Bohr radii
  double step = 0;               // grid spacing, Bohr radii
  std::vector<double> u;         // normalised u(r_i), r_i = i * step
};

// reducedMass is in electron masses. Supports l <= 2.
RadialState solveHydrogenic(int n, int l, double Z, double reducedMass);

// Closed-form non-relativistic energy, -mu Z^2 Ha / (2 n^2), in eV.
double schrodingerEnergyEV(int n, double Z, double reducedMass);

// Dirac-Coulomb energy (point nucleus, no QED) minus rest energy, in eV.
double diracEnergyEV(int n, double j, double Z, double reducedMassMeV);
