# The Limits of Inference in Complex Systems: When Stochastic Models Become Indistinguishable
Fortran implementation of bridge change-of-measure (BCM) and path-integral (PI) techniques for inference in SDEs. 

These codes were used to generate the plots in Ref. [1].

**A Python implementation of these methods is under construction and will be referenced here.**

## Models abbreviations
* OU: Ornstein-Uhlenbeck process (Eq.(1) in Ref.[1]).
* GGM: Generalized Gamma Model (Eq.(2) in Ref.[1]).
* DE: DEmographic model (Eq.(2) in Ref.[1] with $\theta=0$).
* EN: ENvironmental model (Eq.(2) in Ref.[1] with $\theta=1$).
* CP: Contact process (Eq.(A1) in Ref.[1]).
  
## Codes
* **M_declarations.f90**: Define public variables $x$, $t$, $D$, $k$ and $\mu$ to be used in the rest of modules. $x$, and $t$ are used in functions and subroutines to characterize the process at time $t$. $D$, $k$ and $\mu$ are parameters of the models (see Eqs. (1) and (2) in Ref. [1]).
* **M_subroutines.f90**: Main script containing codes for BCM and PI inference methods.
* **M_functions.f90**: Mathematical functions used in the main codes.
* **dranxor.f90**: Pseudo random number generator.
* **GH_draw_trajectory_on_file.f90**: Draw trajectory for selected model and save it to a file (Fig. 4).
* **GH_estimator_of_propagator_BCM.f90**: Estimate propagator using BCM in one transition (Fig.4)
* **GH_model_distinguishability_EN_DE_bridge_CM**: Generate synthetic data and compute model probability using BCM method (Fig. 5).
* **GH_model_distinguishability_EN_DE_PI.f90**: Generate synthetic data and compute model probability using BCM method (Fig. 5).

## Data
* **GH_twizzers.dat** Optical Twizzers data (see License for use). 

## License
This project is shared for **academic and research purposes**. 

The codes are free to use, redistribute, modify, and share for research purposes, provided that proper credit is given to the authors through citation of [1].

Data rights belong to the parties responsible for data acquisition. Any use or distribution of the data should be accompanied by citation of the original data sources (see below) and discussed with the parties responsible for data acquisition.

**Optical twizzers data:** Miguel Ibañez García and Raúl A. Rica, University of Granada. Ref. [2]

## Citation

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

## Reference

[1] Aguilar J., Muñoz M.A., Azaele S. Phys. Rev. X (2026). 10.1103/tmkr-9kl2

[2] Ibáñez, M., Dieball, C., Lasanta, A. et al. Nat. Phys. 20, 135–141 (2024). https://doi.org/10.1038/s41567-023-02269-z
