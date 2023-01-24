___EfficientFrontier.jl___

[![Build Status](https://github.com/PharosAbad/EfficientFrontier.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/PharosAbad/EfficientFrontier.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://github.com/PharosAbad/EfficientFrontier.jl/wiki)
<!-- 
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://PharosAbad.github.io/EfficientFrontier.jl/dev/)
[![Coverage](https://codecov.io/gh/PharosAbad/EfficientFrontier.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/PharosAbad/EfficientFrontier.jl)
-->

<h1 align="center" margin=0px>
  Entire Efficient Frontier by Status-Segment Method
</h1>

<p align="center">
  <a href="#features">Features</a> •
  <a href="#installation">Installation</a> •
  <a href="#license-">License</a> •
  <a href="https://github.com/PharosAbad/EfficientFrontier.jl/wiki">Documentation</a>
</p>

**EfficientFrontier.jl** solves the following problem:

$$
\begin{array}
[c]{cl}
\min & \frac{1}{2}\mathbf{z}^{\prime}\mathbf{Vz}-L\mathbf{z}^{\prime
}\boldsymbol{\mu}\\
s.t. & \mathbf{Az}=\mathbf{b}\in\mathbb{R}^{M}\\
& \mathbf{Gz}\leq\mathbf{g}\in\mathbb{R}^{J}\\
& \boldsymbol{d}\leq\mathbf{z}\leq\boldsymbol{u}\in\mathbb{R}^{N}
\end{array}
$$

with mean vector $\boldsymbol{\mu}\in\mathbb{R}^{N}$ and variance matrix $\mathbf{V}\in\mathbb{R}^{N\times N}$. $\textcolor{blue}{ Varying\ L\ from\ +\infty\ to\ 0}$, all the efficient critical line segments are computed using *closed-form formulas*. And the full Efficient Frontier is recorded by corner portfolios connected by parabolas with *analytical* parameter.

I can't wait to see: [__A Quick Example__](https://github.com/PharosAbad/EfficientFrontier.jl/wiki/A-Quick-Example)

The _Status-Segment Method_ is a two-stage method:
1. find out the _status_ of each asset, that is, whether the weight of the asset falls on the upper or lower bounds of the interval (`OUT`: `DN` or `UP`), or in the middle of the interval (`IN`).  
2. find out the efficient _segment_ of the CL (critical line), that is, the value range of the slope $L$ in the EV plane. 

Since the end points of the efficient segment of a CL provide the status information for the adjacent CLs, the efficient segment of an adjacent CL is found immediately. Therefore, as long as the status of any point on the efficient frontier is found, the entire efficient frontier can be found (one and all).

## Features

* __Entire efficient frontier__: calculate the entire efficient frontier, not just a single frontier portfolio
* __Analytical solutions__: use analytical solutions for calculations, not a numerical method that iterate to convergence (A working paper is coming soon). Even for a single frontier portfolio, our solution is faster and more accurate than the latest numerical methods such as [Clarabel](https://github.com/oxfordcontrol/Clarabel.jl), because we calculate through the analytical solution, rather than through numerical iteration to convergence.
* __Versatile__: from the simplest [no short-sale](https://github.com/PharosAbad/EfficientFrontier.jl/blob/main/examples/frontier.jl) to most general model with [lower and upper bounds, inequality constraints, and equality constraints](https://github.com/PharosAbad/EfficientFrontier.jl/blob/main/examples/ungil.jl). Theoretically we require the variance matrix to be symmetric and positive definite, but in fact we only need the variance matrix subblocks of the IN set to be symmetric and positive definite. Please refer to [SP500](https://github.com/PharosAbad/EfficientFrontier.jl/blob/main/examples/SP500.jl) for a rank-deficient example
* __All-weather__: The [Critical Line Algorithm (CLA)](https://books.google.ch/books?id=eJ8QUsgfZ8wC) (Markowitz, 1956; [Markowitz and Todd, 2000](https://github.com/PharosAbad/EfficientFrontier.jl/blob/main/examples/MarkowitzTodd2000.jl)) requires the model to be non-degenerate (there is only one asset toggling IN and OUT state). Or the perturbation method is used to solve the degenerated cases. Our calculations do not suffer from these problems, and we find [an example](https://github.com/PharosAbad/EfficientFrontier.jl/blob/main/examples/failCLA.jl) of incorrect results obtained by CLA's perturbation algorithm
 * __Plugin__: Simplex method and Combinatorial search method are built-in methods to identify the Status for first CL. [An example](https://github.com/PharosAbad/EfficientFrontier.jl/blob/main/examples/initClarabel.jl) of [plugin](https://github.com/PharosAbad/EfficientFrontier.jl/blob/main/examples/uClarabel.jl) using [Clarabel](https://github.com/oxfordcontrol/Clarabel.jl) is provided, which finds the Status with the highest $L$ value to start the calculation.
* __Open Source__: Our code is available on [GitHub](https://github.com/PharosAbad/EfficientFrontier.jl) and distributed under the MIT License
* __Arbitrary Precision Arithmetic__: fully support for `BigFloat` from v0.2.0

## Installation
__EfficientFrontier.jl__ can be added by

- `import Pkg; Pkg.add("EfficientFrontier")`
- `pkg> add EfficientFrontier`
- `pkg> dev EfficientFrontier` for testing nightly version. To use the registered version again `pkg> free EfficientFrontier`

## License 🔍
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
