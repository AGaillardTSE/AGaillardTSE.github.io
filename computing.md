---
layout: page
---

<h5 style="margin:0 0 10px;padding-left:10px;">
  <a style="color: #000" href="./index.html">Home</a>&nbsp;&nbsp;&nbsp;&nbsp;<a style="color: #000" href="./research.html">Research</a>&nbsp;&nbsp;&nbsp;&nbsp;<a style="color: #000" href="./computing.html">Computing</a>&nbsp;&nbsp;&nbsp;&nbsp;<a style="color: #000" href="./teaching.html">Teaching</a>&nbsp;&nbsp;&nbsp;&nbsp;<a style="color: #000" href="./CV.pdf">Vitae</a>&nbsp;&nbsp;&nbsp;&nbsp;
</h5>
<hr style="margin:0 0 25px;">


# Computing

Below I report some of my codes and replications of recent and former models in macroeconomics, relative to inequality, housing and occupational choice. 

```c++
int C++ // (for discrete time)
```
* Solve standard housing macroeconomic models with DC-EGM algorithm. Useful note here, code available here.
* Solve the Aiyagari model in 0.04 – 0.14 seconds with Endogenous Grid Method (EGM) (Caroll (2006)). Useful note (by Josep Pijoan-Mas) is available here. Download my code (iterate on marginal utilities or value functions with code here).
* Solve the stochastic growth model as in Barillas & Villaverde (2007) using EGM, code here.
* Discretize income process using Tauchen algorithm in C++: code here (by Sumudu Kankanamge

```matlab
MatLab % (for continuous time)
```
* Solve Aiyagari in 0.13 seconds with Envelope Condition Method (ECM), many codes available here: HATC project
* Aiyagari in Continous Time with Jump-Drift Process. Code is available here: aiyagari.m
* Heterogenous Agent New Keynesian (HANK) model and the code available here: (not yet available)

 
### Computational speed - comparison between EGM, DC-EGM and VFI

| Method        | Speed (in s)         | % Entrepreneurs | K/Y |
|:-------------|:------------------|:------|:------|
| EGM           | 0.8s | 8.4  | 2.6 |
| DC-EGM | 1.2s   | 8.8  | 2.6 |
| VFI           | 3s      | 8.8   | 2.6 |



