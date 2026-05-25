# The Limits of Inference in Complex Systems: When Stochastic Models Become Indistinguishable
Codes and data The Limits of Inference in Complex Systems: \\     When Stochastic Models Become Indistinguishable

## Code structure
* **M_declarations.f90**: Define public variables $x$, $y$, $t$, $D$, $r$. $x$, $y$, and $t$ are used in functions and subroutines to characterize the state of the system. $D$ and $r$ represent, respectively, the parameters $D$ and $\rho$ in our text (see License for reference).
* **M_subroutines.f90**: Subroutines like the Milstein method to sample realizations of the process.
* **M_functions.f90**: Mathematical functions used in the main code.
* **dranxor.f90**: Pseudo random number generator.
* **main_Draw_trajectories.f90**: example of code to generate and save $x$, $y$, and $\theta_W$ for different values of $t$.

## License
This project is shared for **academic and research purposes**. It is free to use, redistribute, modify, and share for research purposes, provided that proper credit is given to the authors through citation of: 

Aguilar J., Muñoz M.A., Azaele S. Phys. Rev. X (2026). 10.1103/tmkr-9kl2
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
