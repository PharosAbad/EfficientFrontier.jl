___EfficientFrontier.jl___

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://PharosAbad.github.io/EfficientFrontier.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://PharosAbad.github.io/EfficientFrontier.jl/dev/)
[![Build Status](https://travis-ci.com/PharosAbad/EfficientFrontier.jl.svg?branch=main)](https://travis-ci.com/PharosAbad/EfficientFrontier.jl)
[![Coverage](https://codecov.io/gh/PharosAbad/EfficientFrontier.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/PharosAbad/EfficientFrontier.jl)

<h1 align="center" margin=0px>
  Full Efficient Frontier by connecting Critical Line Segments
</h1>

<p align="center">
  <a href="#features">Features</a> •
  <a href="#installation">Installation</a> •
  <a href="#license-">License</a> •
  <a href="https://PharosAbad.github.io/EfficientFrontier.jl/stable/">Documentation</a>
</p>

EfficientFrontier.jl solves the following problem:

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

with mean vector $\boldsymbol{\mu}\in\mathbb{R}^{N}$ and variance matrix $\mathbf{V}\in\mathbb{R}^{N\times N}$. Varying $L$ from $+\infty$ to $0$, all the efficient ciritlcal line segments are computed using *closed-form formulas*. And the full Efficient Frontier is recorded by corner portfolios connected by parabolas with *analytical* parameter.


## Features

* __Entire efficient frontier__: calculate the entire efficient frontier, not just a single frontier portfolio
* __Analytical solutions__: use analytical solutions for calculations, not a numerical method that iterate to convergence (A working paper is coming soon)
* __Versatile__: from the simplest [no short-sale](blob/main/examples/frontier.jl) to most general model with [lower and upper bounds, inequality constraints, and equality constraints](blob/main/examples/ungil.jl). Theoretically we require the variance matrix to be symmetric and positive definite, but in fact we only need the variance matrix subblocks of the IN set to be symmetric and positive definite. Please refer to [SP500](blob/main/examples/SP500.jl) for a rank-deficient example
* __All-weather__: The CLA algorithm (Markowitz, 1956; Markowitz and Todd, 2000) requires the model to be non-degenerate (there is only one asset toggling IN and OUT state). Or the perturbation method is used to solve the degradation cases. Our calculations do not suffer from these problems, and we find [an example](blob/main/examples/failCLA.jl) of incorrect results obtained by CLA's perturbation algorithm
 * __Novel__: use [Clarabel](https://github.com/oxfordcontrol/Clarabel.jl) to find the critical line with the highest $L$ value to start the calculation. Not only because Clarabel is efficient and fast, but also because the output of Clarabel can provide the position of the weight of each asset in the optimal portfolio in the feasible region, that is, the inner point or the boundary. Which facilitates judging IN, DN, and UP asset collections to initialize critical lines. Conversely, after starting the calculation by Clarabel, our solution to any single frontier combination is faster and more accurate than the best numerical methods such as Clarabel, because we calculate through the analytical solution, rather than through numerical iteration to convergence.
* __Open Source__: Our code is available on [GitHub](https://github.com/PharosAbad/EfficientFrontier.jl) and distributed under the MIT License

## Installation
- __EfficientFrontier.jl__ can be added by `pkg> add https://github.com/PharosAbad/EfficientFrontier.jl.git`

## License 🔍
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
