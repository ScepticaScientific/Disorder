# Disorder/2D
This folder contains files implementing a 2D nonlinear diffusion-reaction model. To start, assemble the project and run `Test.exe`.

### 2D nonlinear diffusion-reaction model
A detailed description of the model, with theoretical results and many numerical experiments, can be found in [1, 2].

### Technical issues
The model was developed in Fortran 90 and tested on Compaq Visual Fortran 6.6 and PGI Visual Fortran 2005. The file `Diffusion.dsp` is the Compaq Visual Fortran project file. If you use this compiler, please make sure you have the IMSL libraries `SMATHD.LIB` and `SF90MP.LIB`.

### References
1. Yu.N. Skiba and D.M. Filatov, On an Efficient Splitting-Based Method for Solving the Diffusion Equation on a Sphere, *Numer. Meth. Part. Diff. Eq.*, 28 (2012) 331-352. <p><a href = "https://doi.org/10.1002/num.20622" rel = "nofollow"><img src = "https://zenodo.org/badge/DOI/10.1002/num.20622.svg" alt = "DOI:10.1002/num.20622" style = "vertical-align: top; max-width: 100%;"></a></p>
2. Yu.N. Skiba and D.M. Filatov, Splitting-based Schemes for Numerical Solution of Nonlinear Diffusion Equations on a Sphere, *Appl. Math. Comput.*, 219 (2013) 8467-8485. <p><a href = "https://doi.org/10.1016/j.amc.2013.02.066" rel = "nofollow"><img src = "https://zenodo.org/badge/DOI/10.1016/j.amc.2013.02.066.svg" alt = "DOI:10.1016/j.amc.2013.02.066" style = "vertical-align: top; max-width: 100%;"></a></p>