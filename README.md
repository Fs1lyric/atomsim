# atomsim

A simulation of an atom, the splitting of a uranium nucleus, and a nuclear chain reaction,
written in C++ and built on real measured physics data.

![Hydrogen orbital](docs/atom.png)
![Uranium fission](docs/fission.png)
![Chain reaction](docs/chain.png)

- **Atom:** the electron cloud of hydrogen, drawn from a solution of the Schrödinger equation.
- **Fission:** a neutron splits uranium-235, and you watch the pieces fly apart.
- **Chain reaction:** neutrons multiply through a ball of uranium. Below the critical mass it
  dies out; above it, it runs away.

The results match real experiments to within a few percent, and the tests check that.

## Run it

```sh
cmake -S . -B build -DATOMSIM_VIZ=ON && cmake --build build -j
./build/atomsim-viz
```

Press 1, 2 or 3 to switch scenes.
