---
layout: page
---


# Computing

Below I report some of my codes and replications of recent and former models in macroeconomics, relative to inequality, housing and occupational choice. 

```c++
int C++ // (for discrete time)
```
* Solve standard housing macroeconomic models (Sommer & Sullivan AER (2018)) with DC-EGM algorithm (Iskhakov et al. (2017)). Useful note <a href="http://agaillard.eu/projects/HOUSING_notes/numerical_solution_Sommer_Sullivan_AER.pdf" target="_blank">here</a>, code available <a href="https://github.com/AGaillardTSE/housing" target="_blank">here</a>.
* Solve the Aiyagari model in 0.04 – 0.14 seconds with Endogenous Grid Method (EGM) (Caroll (2006)). Useful note (by Josep Pijoan-Mas) is available <a href="https://www.cemfi.es/~pijoan/Teaching_files/Notes%20on%20endogenous%20grid%20method.pdf" target="_blank">here</a>. Download my code (iterate on marginal utilities or value functions with code [here](https://github.com/AGaillardTSE/aiyagari)).
* Solve the stochastic growth model as in Barillas & Villaverde (2007) using EGM, code [here](https://github.com/AGaillardTSE/stochastic-growth-model).
* Discretize income process using Tauchen algorithm in C++: code [here].

```matlab
MatLab % (for continuous time)
```
* Solve Aiyagari in 0.13 seconds with Envelope Condition Method (ECM), many codes available here: [HATC project](http://www.princeton.edu/%7Emoll/HACTproject.htm).
* Aiyagari in Continous Time with Jump-Drift Process. Code is available here: aiyagari.m, note: <a href="http://agaillard.eu/resources/aiyagari2.pdf" target="_blank">here</a>.
* Heterogenous Agent New Keynesian (HANK) model and the code available here: (not yet available), note: <a href="http://agaillard.eu/resources/HANK.pdf" target="_blank">here</a>.

 
### Computational speed - comparison between EGM, DC-EGM and VFI
Comparison of performance and accuracy of EGM, DC-EGM and VFI methods on occupational choice and entrepreneurship models à la Cagetti & De Nardi (2006). The presence of discrete choice (occupational choice) makes EGM inaccurate. DC-EGM encompasses generated kinks very well, while being substantially faster than standard VFI.

| Method        | Speed (in s)         | % Entrepreneurs | K/Y |
|:-------------|:------------------|:------|:------|
| EGM           | 0.8s | 8.4  | 2.6 |
| DC-EGM | 1.2s   | 8.8  | 2.6 |
| VFI           | 3s      | 8.8   | 2.6 |



