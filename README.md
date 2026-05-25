# The Limits of Inference in Complex Systems: When Stochastic Models Become Indistinguishable
Fortran implementation of bridge change-of-measure (BCM) and path-integral (PI) techniques for inference in SDEs. These codes were used to generate the plots in Ref. [1].

**A Python implementation of these methods is under construction and will be referenced here.**

## Code structure
* **M_declarations.f90**: Define public variables $x$, $t$, $D$, $k$ and $\mu$ to be used in the rest of modules. $x$, and $t$ are used in functions and subroutines to characterize the process at time $t$. $D$, $k$ and $\mu$ are parameters of the models (see Eqs. (1) and (2) in Ref. [1]).
* **M_subroutines.f90**: Main script containing codes for BCM and PI inference methods.
* **M_functions.f90**: Mathematical functions used in the main codes.
* **dranxor.f90**: Pseudo random number generator.

## License
This project is shared for **academic and research purposes**. It is free to use, redistribute, modify, and share for research purposes, provided that proper credit is given to the authors through citation of [1].

## Reference

[1] Aguilar J., Muñoz M.A., Azaele S. Phys. Rev. X (2026). 10.1103/tmkr-9kl2
```
@article{tmkr-9kl2,
  title = {Limits of inference in complex systems: When stochastic models become indistinguishable},
  author = {Aguilar, Javier and Munoz, Miguel A. and Azaele, Sandro},
  journal = {Phys. Rev. X},
  pages = {},
  year = {2026},
  month = {May},
  publisher = {American Physical Society},
  doi = {10.1103/tmkr-9kl2},
  url = {https://link.aps.org/doi/10.1103/tmkr-9kl2}
}
```
