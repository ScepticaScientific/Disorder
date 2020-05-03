# Disorder/3D
This folder contains files implementing a 3D nonlinear diffusion-reaction model. To start, assemble the project and run `Test.exe`.

### 3D nonlinear diffusion-reaction model
A detailed description of the model, with theoretical results and many numerical experiments, can be found in [1-2].

### Technical issues
The Python Fourier CCA code involves calling the function `getPowerSpectrum()` from `../utils`. The code of `getPowerSpectrum()` ensures an accurate computation of power spectra of multivariate datasets regardless of the length of the time series. The function `getPowerSpectrum()` was developed once the standard Python function `scipy.signal.csd()` had been found to provide inaccurate results on long time series containing 10,000+ samples each.

### References
1. Yu.N. Skiba and D.M. Filatov, Modelling of Combustion and Diverse Blow-Up Regimes in a Spherical Shell, in: P. Quintela *et al.* (eds.), *Rev. Sel. Papers of the 19th European Conf. on Mathematics for Industry (ECMI 2016), Progress in Industrial Mathematics at ECMI 2016, Mathematics in Industry*, Springer, Cham, 2017, pp. 729-735. <p><a href = "https://doi.org/10.1007/978-3-319-63082-3_108" rel = "nofollow"><img src = "https://zenodo.org/badge/DOI/10.1007/978-3-319-63082-3_108.svg" alt = "DOI:10.1007/978-3-319-63082-3_108" style = "vertical-align: top; max-width: 100%;"></a></p>
2. Yu.N. Skiba and D.M. Filatov, A Numerical Study of Nonlinear Diffusion Phenomena in Heterogeneous Media: Energy Transfer at Diverse Blow-Up Modes and Self-Organisation Processes, *Eur. Phys. J. Special Topics*, 226 (2017) 3303-3314. <p><a href = "https://doi.org/10.1140/epjst/e2016-60323-x" rel = "nofollow"><img src = "https://zenodo.org/badge/DOI/10.1140/epjst/e2016-60323-x.svg" alt = "DOI:10.1140/epjst/e2016-60323-x" style = "vertical-align: top; max-width: 100%;"></a></p>
