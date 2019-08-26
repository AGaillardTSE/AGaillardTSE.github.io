---
layout: page
---


# Computing

##### Useful papers & technics:
* Endogenous Grid Method (EGM) ([Caroll (2006)](http://pages.stern.nyu.edu/~dbackus/Computation/Carroll%20endog%20grid%20EL%2006.pdf))
* DC-EGM algorithm to solve discrete-continuous choice models using EGM ([Iskhakov et al. (2017)](https://onlinelibrary.wiley.com/doi/abs/10.3982/QE643)).
* Computing expectations of value functions using polynomials: [Judd et al. (2017)](https://onlinelibrary.wiley.com/doi/abs/10.3982/QE329) 
* Multi-dimensional DC-EGM algorithm: [Druedhal and Jørgensen (2017)](https://www.sciencedirect.com/science/article/pii/S0165188916301920)

Below I report some of my codes and replications of recent and former models in macroeconomics, relative to inequality, housing and occupational choice. 

```c++
int C++ // (for discrete time)
```
* Solve standard housing macroeconomic models ([Sommer & Sullivan AER (2018)](https://www.aeaweb.org/articles?id=10.1257/aer.20141751)) with DC-EGM. Useful note <a href="http://agaillard.eu/projects/HOUSING_notes/numerical_solution_Sommer_Sullivan_AER.pdf" target="_blank">here</a>, code available <a href="https://github.com/AGaillardTSE/housing" target="_blank">here</a>.
* Solve the Aiyagari model in 0.04 – 0.14 seconds with EGM. Useful note (by Josep Pijoan-Mas) is available <a href="https://www.cemfi.es/~pijoan/Teaching_files/Notes%20on%20endogenous%20grid%20method.pdf" target="_blank">here</a>. Download my code (iterate on marginal utilities or value functions with code [here](https://github.com/AGaillardTSE/aiyagari)).
* Solve the stochastic growth model as in [Barillas & Villaverde (2007)](https://econpapers.repec.org/article/eeedyncon/v_3a31_3ay_3a2007_3ai_3a8_3ap_3a2698-2712.htm) using EGM, code [here](https://github.com/AGaillardTSE/stochastic-growth-model).
* Discretize income process using Tauchen algorithm in C++: code [here].

```matlab
MatLab % (for continuous time)
```
* Solve Aiyagari in 0.13 seconds with Envelope Condition Method (ECM), many codes available here: [HATC project](http://www.princeton.edu/%7Emoll/HACTproject.htm).
* Aiyagari in Continous Time with Jump-Drift Process. Code is available here: aiyagari.m, note: <a href="http://agaillard.eu/resources/aiyagari2.pdf" target="_blank">here</a>.
* Heterogenous Agent New Keynesian (HANK) ([Kaplan et al. (2018)](https://www.princeton.edu/~moll/HANK.pdf)) model and the code available here: (not yet available), note: <a href="http://agaillard.eu/resources/HANK.pdf" target="_blank">here</a>.

 
### Computational speed - comparison between EGM, DC-EGM and VFI
Comparison of performance and accuracy of EGM, DC-EGM and VFI methods on occupational choice and entrepreneurship models à la [Cagetti & De Nardi (2006)](http://users.nber.org/~denardim/research/caciocristina.pdf). The presence of discrete choice (occupational choice) makes EGM inaccurate. DC-EGM encompasses generated kinks very well, while being substantially faster than standard VFI.

| Method        | Speed (in s)         | % Entrepreneurs | K/Y |
|:-------------|:------------------|:------|:------|
| EGM           | 0.8s | 8.4  | 2.6 |
| DC-EGM | 1.2s   | 8.8  | 2.6 |
| VFI           | 3s      | 8.8   | 2.6 |




## Useful Links ##
* [Jean-Pierre's Moreau Homepage](http://jean-pierre.moreau.pagesperso-orange.fr/links.html): useful codes and routines in C++ and Fortran.
* [John Starchulski and Tom Sargent's QuantEcon](http://quant-econ.net): useful codes in Julia / Python.




