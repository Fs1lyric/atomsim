# atomsim

An atom, the splitting of its nucleus, and a nuclear chain reaction, simulated in C++17
from measured nuclear data and checked against experiment in CI.

| The atom | One fission | Chain reaction |
|---|---|---|
| ![Hydrogen 3d orbital](docs/atom.png) | ![U-235 fission](docs/fission.png) | ![Supercritical uranium sphere](docs/chain.png) |
| Electron cloud sampled from a numerical solution of the Schrödinger equation | ²³⁶U* splits; fragments accelerate along the computed Coulomb trajectory, then emit neutrons | Neutron transport through uranium with ENDF/B-VIII.0 cross-sections |

## Results vs experiment

| Quantity | Model | Measured | Off by |
|---|---|---|---|
| Godiva critical mass (bare 93.7% U-235 sphere) | 52.92 kg | 52.42 kg | 1.0% |
| Godiva k_eff at its measured size | 0.9993 | 1.0000 | 0.07% |
| Prompt neutrons per U-235 fission ν | 2.33 | 2.41 | 3.5% |
| Fragment kinetic energy | 171.7 MeV | 170.5 MeV | 0.7% |
| Light / heavy fragment mass | 96.8 / 139.2 | 96.4 / 139.6 | <0.5% |
| Prompt neutron energy | 1.92 MeV | 2.0 MeV | 4% |
| Prompt γ energy | 6.74 MeV | 6.6–6.9 MeV | in range |
| Hydrogen 1s energy (Numerov vs exact) | −13.598287264 eV | same | 2×10⁻¹³ eV |

Every row is a test in `tests/test_physics.cpp` (GoogleTest). CI fails if any of them
drifts outside its stated tolerance.

## Build and run

```sh
cmake -S . -B build -DATOMSIM_VIZ=ON && cmake --build build -j
./build/atomsim-viz                      # the visualisation (keys 1-3 switch scenes)
./build/atomsim                          # full text report: atom, fission, Monte Carlo
./build/atomsim --target U238 --energy 2 # fast-neutron fission of U-238
ctest --test-dir build                   # physics tests
```

The visualisation also builds to WebAssembly with Emscripten (see the `web` CI job).

## How it works

```
src/hydrogen   Numerov integration + shooting for the radial Schrödinger equation; Dirac energies
src/nucleus    AME2020 masses (3,500 nuclides), Myers–Swiatecki liquid drop with shell corrections
               as fallback, fission barrier, Hill–Wheeler transmission
src/fission    one fission event: fragment sampling, Q-value, Coulomb acceleration (relativistic
               RK4), fragment spin, neutron evaporation, γ cascade, β chains; parallel Monte Carlo
src/chain      continuous-energy neutron transport: elastic scattering with ENDF angular
               distributions, inelastic levels, (n,2n), fission, capture; k_eff and time evolution
viz/           raylib front end; every number on screen comes from the library above
tools/         extract_endf turns 50 MB of ENDF-6 evaluations into 230 KB of tables
```

## Design decisions

- **Measured data over fitted formulas.** Fragment masses come from the AME2020 table, not
  the textbook mass formula. The formula was off by up to 10 MeV near magic numbers and
  pushed ν to 3.2; measured masses brought it to 2.33.
- **Continuous-energy transport over one-group.** A one-group model with textbook constants
  predicted a 38 kg critical mass (27% low). Neutrons scattering off uranium at MeV energies
  mostly go forwards, and inelastic collisions slow them below U-238's fission threshold.
  Using the ENDF angular distributions and inelastic levels fixed it (1% off).
- **Fragment spin.** Without it, ν was too high and γ energy too low. Fragments spin at about
  7ħ, and that rotational energy can only leave as γ rays. Adding it fixed both at once.
- **Deterministic parallelism.** The Monte Carlo splits work into independently seeded
  blocks merged in a fixed order, so results are bit-identical on 1 or 64 threads.
  A test enforces this.
- **Rendering is honest about its limits.** The nucleus stretching before scission is drawn
  schematically and the app says so; everything after scission is computed.

## Limits

- The fission mass distribution (Brosa modes) and the effective scission distance are
  fitted to data, not derived. Deriving them needs Langevin dynamics on a computed
  potential-energy surface; `ROADMAP.md` describes how.
- Below 25 keV, transport uses the 25 keV cross-sections. Resolved resonances are not
  reconstructed, which matters little for fast systems like Godiva.
- Not modelled: QED, delayed neutrons, half-lives, thermal-neutron reactors.

## Data sources

- Nuclear masses: AME2020, Wang et al., *Chinese Phys. C* 45, 030003 (2021)
- Cross-sections: ENDF/B-VIII.0, Brown et al., *Nucl. Data Sheets* 148, 1 (2018)
- Benchmark: ICSBEP HEU-MET-FAST-001 (Godiva)
- Fission observables: Madland, *Nucl. Phys. A* 772 (2006); ENDF/B-VIII.0 ν̄
