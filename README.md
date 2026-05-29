# MPS/MPO Chebyshev Expansion — Density of States via Tensor Networks

A MATLAB implementation of the **Kernel Polynomial Method (KPM)** using **Matrix Product States (MPS)** and **Matrix Product Operators (MPO)** to compute the density of states (DOS) of quantum lattice Hamiltonians.

The primary application is the transverse-field Ising model, but the framework is general and can be adapted to any 1D Hamiltonian expressible as an MPO.

---

## Repository Structure

```
mps-chebyshev/
├── src/                        # All source functions
│   ├── ChebyshevH.m            # Main: builds Chebyshev MPO expansion of H
│   ├── trChebyshevH.m          # Computes Tr[T_n(H)] or Tr[O * T_n(H)]
│   ├── dos.m                   # Density of states via KPM (standard)
│   ├── dos_precise.m           # DOS with separate spectral resolution
│   ├── reduceD.m               # Variational MPO/MPS bond-dimension compression
│   ├── minimizeE.m             # DMRG ground-state energy minimization
│   ├── maximizeE.m             # DMRG maximum eigenstate search
│   ├── prepare.m               # Right-to-left MPS canonicalization
│   ├── prepare_onesite.m       # Single-site SVD gauge fixing
│   ├── createrandommps.m       # Random MPS initialization
│   ├── initCstorage.m          # Initialize environment tensors (C storage)
│   ├── updateCleft.m           # Left-to-right environment update
│   ├── updateCright.m          # Right-to-left environment update
│   ├── contracttensors.m       # Generic tensor contraction engine
│   ├── HSnorm.m                # Hilbert-Schmidt inner product of two MPOs
│   ├── svd2.m                  # Memory-efficient thin SVD wrapper
│   └── project_H.m             # Spectral projection helper
├── examples/
│   ├── ising.m                 # DMRG ground/max state of Ising model
│   └── isingtr.m               # Full DOS calculation for Ising model
├── docs/
│   └── method_notes.md         # Extended notes on the algorithm
└── README.md
```

---

## Quickstart

### 1. Ground state energy (DMRG)

```matlab
cd examples
ising   % runs DMRG for N=40 Ising chain, D=14 bond dimension
```

Output: `E0` (ground state energy), `E1` (maximum eigenvalue energy).

### 2. Density of states (Chebyshev + KPM)

```matlab
cd examples
isingtr  % computes DOS and <sigma_z> spectral function for N=40
```

Output files: `trH.txt`, `trOH.txt`.  
A plot of $\langle \sigma^z_{N/2} \rangle$ vs energy is displayed automatically.

---


