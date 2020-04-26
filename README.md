## Projected Interleaved Agile Factor Decomposition (PIAFD)

This code accompanies the following publication:

Bai J, Ament S, Perez G, Gregoire J, Gomes C. ***An Efficient Relaxed Projection Method for Constrained Non-negative Matrix Factorization with Application to the Phase-Mapping Problem in Materials Science.*** InInternational Conference on the Integration of Constraint Programming, Artificial Intelligence, and Operations Research 2018 (pp. 52-62). Springer, Cham.

Please cite this work as appropriate in follow-up publications.

Bibtex:

```
@inproceedings{bai2018efficient,
  title={An Efficient Relaxed Projection Method for Constrained Non-negative Matrix Factorization with Application to the Phase-Mapping Problem in Materials Science},
  author={Bai, Junwen and Ament, Sebastian and Perez, Guillaume and Gregoire, John and Gomes, Carla},
  booktitle={International Conference on the Integration of Constraint Programming, Artificial Intelligence, and Operations Research},
  pages={52--62},
  year={2018},
  organization={Springer}
}
```



### Introduction and Prerequisites

PIAFD is an Non-negative Matrix Factorization (NMF) based demixing method. Compared to [AgileFD](https://github.com/JunwenBai/AgileFD), [IAFD](https://github.com/JunwenBai/IAFD), PIAFD gets rid of the need of CPLEX  Optimization Studio and greatly saves running time by relaxing the original phase-mapping problem to a softly constrained optimization problem , which could be optimized by gradient descents. Meanwhile, PIAFD could roughly maintain the same performance as IAFD.

PIAFD has only 1 external library dependencies which must be installed:

1) Armadillo, see: http://arma.sourceforge.net

### Compilation

An example Makefile for GNU make and g++ is included, but will need to be modified to reflect the installation location of the BLAS library (like openblas) used with Armadillo. Additional changes will be needed for alternate compilers. This configuration has been used on a few different Linux distributions, but additional changes could be required for other systems.

After installing all the prerequisites, one should be able to compile the code with the following command:

```bash
make
```

If the compilation fails, try to `make clean` first.

Also, a pre-compiled `solver` is available. One can skip the compilation and move on to usage directly.

### Usage

```bash
./solver input/Ta-Rh-Pd/config.txt

```

Or

```bash
./run.sh
```

One can change the `Ta-Rh-Pd` to any other material system of interest.

### Contact me

If you have any questions or suggestions, please feel free to contact jb2467@cornell.edu or rab38@cornell.edu.
