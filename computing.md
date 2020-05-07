---
layout: page
---


# Computing HA models

#### Useful papers & technics:
##### Endogenous Grid Method (EGM)
 - **[EGM]** Endogenous Grid Method: [Caroll (2006)](http://pages.stern.nyu.edu/~dbackus/Computation/Carroll%20endog%20grid%20EL%2006.pdf)
 - **[DC-EGM]** algorithm to solve discrete-continuous choice models using EGM: [Iskhakov et al. (2017)](https://onlinelibrary.wiley.com/doi/abs/10.3982/QE643). *[G2-EGM]* Multi-dimensional DC-EGM algorithm: [Druedhal and Jørgensen (2017)](https://www.sciencedirect.com/science/article/pii/S0165188916301920). *[N-EGM]* Nested Engodenous Grid Method :: [Hintermaier and Koeniger (2012)](https://hal.archives-ouvertes.fr/hal-00732758/document), [Druedhal (2019)](http://web.econ.ku.dk/druedahl/papers/2019_NEGM.pdf)
 
##### Other technics
 
* **[EXPECTATION 1]** Useful tool to approximate expectations over normally distributed shocks: Gauss-Hermite + Monomial Rule.
* **[EXPECTATION 2]** Computing expectations of value functions using polynomials: [Judd et al. (2017)](https://onlinelibrary.wiley.com/doi/abs/10.3982/QE329)  
* **[SIMULATION]** Non-stochastic simulation routine: [Young (2010)](http://people.virginia.edu/~ey2d/young_2010.pdf)
* **[COLOCATION]** Collocation Method: [Collocation](https://en.wikipedia.org/wiki/Collocation_method)
* **[ECM]** Envelop Condition Method (ECM): [Maliar (2013)](https://stanford.edu/~maliarl/Files/EL2013.pdf)
* **[SMM]** Estimating HA models :: [Sobol sequence](https://en.wikipedia.org/wiki/Sobol_sequence) 
* **[VAR, AR(1)]** Discretizing VAR: [Farmer et Toda (2013)](https://www.econstor.eu/bitstream/10419/195551/1/1015515150.pdf)  

#### Replications

```c++
int C++ // (for discrete time)
```
* **[HOUSING]** Solve standard housing macroeconomic models with DC-EGM, [Sommer & Sullivan AER (2018)](https://www.aeaweb.org/articles?id=10.1257/aer.20141751). Useful note <a href="http://agaillard.eu/projects/HOUSING_notes/numerical_solution_Sommer_Sullivan_AER.pdf" target="_blank">here</a>, code available *upon-request*.
* **[AIYAGARI]** Solve the Aiyagari model in 0.04 – 0.14 seconds with EGM + Young method (2010). Useful note (by Josep Pijoan-Mas) is available <a href="https://www.cemfi.es/~pijoan/Teaching_files/Notes%20on%20endogenous%20grid%20method.pdf" target="_blank">here</a>. Download my code (iterate on marginal utilities or value functions with code available *upon-request*).
* **[SGM]** Solve the stochastic growth model using EGM, [Barillas & Villaverde (2007)](https://econpapers.repec.org/article/eeedyncon/v_3a31_3ay_3a2007_3ai_3a8_3ap_3a2698-2712.htm), code [here](https://github.com/AGaillardTSE/stochastic-growth-model).
* **[ENTREPRENEURSHIP]** Solve [Cagetti & DeNardi (2006)](http://users.nber.org/~denardim/research/JPEfinal.pdf) in 3s with DC-EGM, code available *upon-request*.

```matlab
MatLab % (for continuous time)
```
* **[AIYAGARI]** Solve Aiyagari in 0.13 seconds with Envelope Condition Method (ECM), many codes available here: [HATC project](http://www.princeton.edu/%7Emoll/HACTproject.htm).
* **[AIYAGARI]** Aiyagari in Continous Time with Jump-Drift Process. Code is available *upon-request*, note: <a href="http://agaillard.eu/resources/aiyagari2.pdf" target="_blank">here</a>.
* **[HANK]** Heterogenous Agent New Keynesian (HANK) ([Kaplan et al. (2018)](https://www.princeton.edu/~moll/HANK.pdf)) model and the code available here: (not yet available), note: <a href="http://agaillard.eu/resources/HANK.pdf" target="_blank">here</a>.

 
#### Computational speed - comparison between EGM, DC-EGM and VFI
Comparison of performance and accuracy of EGM, DC-EGM and VFI methods on occupational choice and entrepreneurship models à la [Cagetti & De Nardi (2006)](http://users.nber.org/~denardim/research/caciocristina.pdf). The presence of discrete choice (occupational choice) makes EGM inaccurate. DC-EGM encompasses generated kinks very well, while being substantially faster than standard VFI.

| Method        | Speed (in s)         | % Entrepreneurs | K/Y |
|:-------------|:------------------|:------|:------|
| EGM           | 0.8s | 8.4  | 2.6 |
| DC-EGM | 1.2s   | 8.8  | 2.6 |
| VFI           | 3s      | 8.8   | 2.6 |




#### Useful Links 
* [Jean-Pierre's Moreau Homepage](http://jean-pierre.moreau.pagesperso-orange.fr/links.html): useful codes and routines in C++ and Fortran.
* [John Starchulski and Tom Sargent's QuantEcon](http://quant-econ.net): useful codes in Julia / Python.
* [Paul Bourke](http://paulbourke.org/miscellaneous/interpolation/): useful tool for interpolations.
* [John Burkardt](https://people.sc.fsu.edu/~jburkardt/): useful tool for many languages.
* [Maliar et Maliar's Book](https://web.stanford.edu/~maliarl/Files/Chapter7HandbookCE2014.pdf): all.


