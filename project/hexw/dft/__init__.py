"""hexw.dft: Quantum ESPRESSO adapter - collinear magnetic configurations, exchange constants J of a Heisenberg model.

Theory-side code: the vault and holdout are closed to it by hexw.io.paths, programs start only through hexw.util.external.run_external.
Convention: H = -sum_{i<j} J_ij S_i.S_j (|S| spin quantum numbers, J > 0 ferromagnetic, meV); see docs/DFT.md and hexw.dft.hamiltonian.
Importing the package does not read data or start anything.
"""
