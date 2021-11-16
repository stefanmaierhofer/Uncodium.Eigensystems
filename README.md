[![windows](https://github.com/stefanmaierhofer/Uncodium.Eigensystems/actions/workflows/windows.yml/badge.svg)](https://github.com/stefanmaierhofer/Uncodium.Eigensystems/actions/workflows/windows.yml)
[![linux](https://github.com/stefanmaierhofer/Uncodium.Eigensystems/actions/workflows/linux.yml/badge.svg)](https://github.com/stefanmaierhofer/Uncodium.Eigensystems/actions/workflows/linux.yml)
[![macos](https://github.com/stefanmaierhofer/Uncodium.Eigensystems/actions/workflows/macos.yml/badge.svg)](https://github.com/stefanmaierhofer/Uncodium.Eigensystems/actions/workflowsmacoswindows.yml)

# Efficient computation of eigenvalues and eigenvectors of 3x3 matrices.
See original work at https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/.

A common scientific problem is the numerical calculation of the
eigensystem of symmetric or hermitian 3x3 matrices. If this
calculation has to be performed many times, standard packages
like LAPACK, the GNU Scientific Library, and the Numerical Recipes
Library may not be the optimal choice because they are optimized
mainly for large matrices.

This is a C# port of the C implementations of several algorithms
which were specifically optimized for 3x3 problems.
All algorithms are discussed in detail in the following paper:

```
Joachim Kopp
Efficient numerical diagonalization of hermitian 3x3 matrices
Int.J.Mod.Phys.C 19 (2008) 523-548
arXiv.org: http://arxiv.org/abs/physics/0610206
```