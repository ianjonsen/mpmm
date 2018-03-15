# mpmm
(animal) Move Persistence Mixed-effects Models

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

This is a development package that fits a random walk with time-varying move persistence (autocorrelation) to multiple individual animal tracks. Move persistence (gamma_t) is a latent variable that is estimated as a linear mixed-effects function of user-supplied covariates. The models are fitted using maximum likelihood estimation via the `TMB` ([Template Model Builder](https://github.com/kaskr/adcomp)) package and using C++. Random effects are assumed to be Gaussian. The `mpmm` package leverages tools from the `lme4` and `glmmTMB` packages, adopting their formula syntax for random effects components. 

## Installation 
Currently, for testing and evaluation purposes only. You will need C++ compile tools. 

On PC's running Windows, ensure you have installed [Rtools](https://cran.r-project.org/bin/windows/Rtools/) 

On Mac's, ensure you have installed [Xcode](https://developer.apple.com/xcode/) and Xcode developer tools. If installation is needed, make sure you start Xcode after install to ensure final setup of developer tools is completed. Both Xcode and Xcode developer tools can be installed from the [Mac App Store](https://itunes.apple.com/au/app/xcode/id497799835?mt=12)

### From GitHub
`mpmm` development version is available via:
```
devtools::install_github("ianjonsen/mpmm")
```
