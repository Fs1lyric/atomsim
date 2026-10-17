# atomsim

An atom, and the splitting of its nucleus, simulated in C++17 from measured data and
established theory. No dependencies; output goes to the terminal.

```
cmake -S . -B build && cmake --build build
./build/atomsim                          # U-235 + thermal neutron
./build/atomsim --target U238            # below the fission barrier: no split
./build/atomsim --target U238 --energy 2 # fast neutron: it splits
./build/atomsim --target Pu239 --events 50000 --seed 1
ctest --test-dir build                   # physics checks
```

## What it computes

**Part 1: the atom**
- Hydrogen: the radial Schrödinger equation is integrated with Numerov's method and the
  energies are found by shooting. They match the exact `-μ/2n²` Hartree to ~10⁻⁹ eV.
  Dirac energies give the fine structure (2p₃/₂–2p₁/₂ = 10.97 GHz).
- The heavy atom: its binding energy from AME2020 masses vs the semi-empirical mass
  formula (SEMF), its size, its electron configuration (Madelung rule plus the measured
  exceptions), and the binding-per-nucleon curve that makes fission release energy.

**Part 2: fission**
1. **Capture.** The neutron forms a compound nucleus with excitation `E* = S_n + E_n·A/(A+1)`,
   where S_n comes from AME2020 masses.
2. **Barrier.** The liquid-drop fission barrier comes from the Cohen–Swiatecki fit with
   Myers–Swiatecki fissility, and the tunnelling probability from Hill–Wheeler. This is
   why U-235 splits with slow neutrons and U-238 doesn't: S_n is 6.55 vs 4.81 MeV because
   of neutron pairing.
3. **Split.** Fragment masses are sampled from Brosa's fission modes (Standard I near
   ¹³²Sn, Standard II, and the symmetric Superlong mode). Charges follow the unchanged-charge
   distribution with Wahl charge polarisation.
4. **Energy.** `Q = B(light) + B(heavy) − B(parent)`. TKE is the Coulomb energy at
   scission, and the rest becomes fragment excitation.
5. **Acceleration.** The relativistic two-body Coulomb motion is integrated with RK4 from
   ~10⁻²⁴ s to ~10⁻¹⁷ s.
6. **Prompt neutrons.** Weisskopf evaporation with Fermi-gas temperatures, boosted into
   the lab frame. The remaining energy goes out as prompt γ rays.
7. **β⁻ chains.** Each fragment decays to its stable isobar, with SEMF Q-values.
8. **Monte Carlo.** Mass yields, P(ν), the neutron spectrum (compared with the evaluated
   Watt spectrum) and the energy budget, all compared with measured U-235(n_th,f) data.

## Limits (honest list)
- Fragment masses use the Myers–Swiatecki liquid drop with spherical shell corrections;
  strongly deformed fragments are still off by a few MeV.
- Mass-yield modes and the effective scission distance (1.76 fm·(A_L^⅓+A_H^⅓)) are
  fitted to data, not derived. Deriving them needs time-dependent density-functional
  theory on a supercomputer.
- Fragment spins (mean 7ħ) use a rigid-rotor moment of inertia. Not modelled: QED (Lamb
  shift), the sawtooth in excitation sharing, delayed neutrons,
  and half-lives.
