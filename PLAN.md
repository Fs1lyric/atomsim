# Plan: derive everything that is currently fitted

Goal: replace every tuned input with a calculation from established theory, and close
the remaining gaps with experiment (prompt ν is 2.64 in the model vs 2.41 measured).
"Exact" in the strict sense is impossible: fission has no exact solution, and the best
microscopic methods (time-dependent DFT) need supercomputers. The target here is the
method used by real fission codes (Möller/Randrup, FREYA, CGMF): macroscopic-microscopic
energies, Langevin dynamics and Hauser–Feshbach decay.

Rules: all builds and tests run in GitHub Actions, never locally. Write the test first,
watch it fail in CI, then implement. Each phase ships on its own and keeps CI green.

## Current fitted inputs, and what replaces each one
| Fitted now | Replaced by | Phase |
|---|---|---|
| Spherical Myers–Swiatecki shell term | Strutinsky shell + BCS pairing from a deformed Nilsson spectrum, minimised over (β2, β4) | 1 |
| Cohen–Swiatecki barrier, ħω = 0.6 | Barrier and curvature read off the computed potential-energy surface (double hump) | 2 |
| Brosa mass modes (weights, means, widths) | Langevin trajectories on the surface; yields emerge from where they cut | 3 |
| Scission distance 1.76 fm·(A_L^⅓+A_H^⅓), TKE width 8 MeV | Coulomb + pre-scission kinetic energy of the actual scission shape | 3 |
| Excitation shared by mass | Deformation energy of each fragment at scission + equal temperature (gives the sawtooth) | 3 |
| Charge polarisation 0.5, width 0.56 | Minimise fragment-pair energy in Z at fixed A; thermal width from T | 3 |
| Spin cutoff B = 6 | σ² = 𝓘T/ħ² from bending/wriggling at the scission temperature | 3 |
| Level density a = A/11 | Ignatyuk shell-damped level density from Phase 1 shell energies | 4 |
| "Evaporate until E* < S_n" | Hauser–Feshbach n/γ competition: Koning–Delaroche optical-model transmissions, GLO γ strength | 4 |
| No fission/capture branching | Γ_f (Bohr–Wheeler) vs Γ_n, Γ_γ → σ_f/σ_abs (U-235 thermal ≈ 0.85) | 4 |
| β chains: Q only | Fermi theory spectra (Coulomb-corrected) for e⁻/ν̄ split; half-lives from ENSDF data fetched in CI | 5 |
| Hydrogen: no QED | Leading Lamb shift (Bethe logarithm) | 6 |
| Uranium electrons: configuration only | Self-consistent Dirac–Fock–Slater central field; K-shell binding vs measured 115.6 keV | 6 |

## Phases and how each is verified
1. **Deformed mass model.** `src/nilsson.*` (modified oscillator, Bengtsson–Ragnarsson
   parameters, diagonalised per N-shell), `src/strutinsky.*` (smoothing order 6, γ = 1.2ħω),
   BCS pairing. In CI, fetch AME2020 (`mass.mas20` from the IAEA AMDC). Tests: rms error
   < 1.0 MeV over fission-fragment nuclides (Z 30–64), and deformed ground states for
   ¹⁰⁰Zr and ¹⁵²Nd (β2 ≈ 0.3).
2. **Fission surface.** Axial shapes in (elongation, neck, mass asymmetry). Surface and
   Coulomb energies by numerical integration; shell correction from Phase 1 on the same
   shapes. Compute it once in CI and cache it as `data/pes-U236.bin` (artifact, then
   committed). Tests: ²³⁶U inner and outer barriers within 1 MeV of 5.6 and 5.5 MeV, and
   an asymmetric saddle.
3. **Langevin dynamics.** 3D Langevin with Werner–Wheeler inertia, wall-and-window
   friction, and fluctuation–dissipation noise at T(E*). Scission when the neck radius
   is under about 1.5 fm. Tests: ⟨A_H⟩ = 139.6 ± 1, peak width within 15%, ⟨TKE⟩ within
   3 MeV of 170.5, and the sawtooth ν(A) has the right sign. Delete `kStandard1/2`,
   `kScissionRadius`, `kTkeWidth`, `kPolarisation`, `kChargeWidth`, `kSpinCutoff`.
4. **Statistical decay.** Hauser–Feshbach cascade per fragment. Tests: ν = 2.41 ± 0.1,
   P(ν) within 0.03 per bin, γ energy 6.6–6.9, U-235 σ_f/σ_abs ≈ 0.85.
5. **β decay.** Tests: per-fission β energy ≈ 7 MeV, ν̄ energy ≈ 10 MeV.
6. **Atomic physics.** Tests: hydrogen 1S Lamb shift ≈ 8.17 GHz; U K-shell within 1%.

Runtime budget: the surface is the expensive step (minutes in CI, cached). Monte Carlo
stays under a minute.

## Status
- [x] Current version: all checks green (commit b71ee6a)
- [ ] Phase 1 — start here
