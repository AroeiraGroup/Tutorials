## Introduction

In the [first tutorial](../1_MonteCarlo/instructions.md), where we estimated the value of $\pi$, you were introduced to the basic idea of Monte Carlo methods through simple random sampling. In particular, you explored how geometric quantities such as $\pi$ can be estimated by interpreting them as expectations over random variables.

*Monte Carlo Integration* generalizes this idea. Instead of estimating a geometric probability, we aim to compute definite integrals by expressing them as statistical averages. For an integral of the form

$$
I = \int_D f(x)\, dx
$$

*In its simplest form*, Monte Carlo integration approximates $I$ by sampling points $x_i$ uniformly from the domain $D$ and computing

$$
I \approx \frac{|D|}{N} \sum_{i=1}^N f(x_i)
$$

where $|D|$ is the volume of the integration domain and $N$ is the number of samples.

Just as in the estimation of $\pi$, the result of a Monte Carlo integration is a random variable. Its statistical error decreases as $1/\sqrt{N}$, independently of the dimension of the problem. This dimension-independent convergence makes Monte Carlo methods particularly useful for high-dimensional integrals, where deterministic quadrature methods quickly become inefficient.

In this tutorial, you will build on the concepts introduced in the $pi$ estimator problem by applying Monte Carlo sampling to numerical integration. You will start with plain Monte Carlo estimators in one dimension and then introduce variance-reduction techniques that dramatically improve accuracy. Finally, you will explore extensions to higher-dimensional integrals and discuss the strengths and limitations of the method.

## First task: evaluate a one-dimensional integral

A good first exercise is to estimate an integral that is nontrivial but has a simple, exact value so students can check correctness and study error statistics. We use

$$
I = \int_0^{\pi} x \sin x \, dx.
$$

This integral can be evaluated analytically by integration by parts:

$$
\int_0^{\pi} x\sin x\,dx = \big[-x\cos x + \sin x\big]_0^{\pi} = \pi.
$$

So the exact value is $I=\pi$.

### Monte Carlo estimator (uniform sampling)

In this first exercise we use *uniform sampling* on the interval $[0,\pi]$. That is, the random variables $X_i$ are drawn from the probability density

$$ p(x) = \frac{1}{\pi}, \qquad x \in [0,\pi]. $$

With this choice, the integral can be written as

$$
I = \int_0^{\pi} x\sin x\,dx
= \int_0^{\pi} \frac{x\sin x}{p(x)}\,p(x)\,dx
= \mathbb{E}_p\!\left[\pi\,x\sin x\right].
$$

### Error bars

Let $X_1,\dots,X_N$ be independent samples drawn uniformly from $[0,\pi]$, i.e.
$p(x)=1/\pi$, and define the Monte Carlo weights
$$
w_i = \pi\, f(X_i) = \pi\, X_i \sin X_i.
$$

The Monte Carlo estimator of the integral is

$$
\hat I_N = \frac{1}{N}\sum_{i=1}^N w_i.
$$

An estimate of the statistical uncertainty (standard error) of $\hat I_N$ is

$$
\mathrm{SE}(\hat I_N)
=
\frac{s}{\sqrt{N}},
$$

where

$$
s^2
=
\frac{1}{N-1}
\sum_{i=1}^N
\left(w_i - \bar w\right)^2,
\qquad
\bar w = \frac{1}{N}\sum_{i=1}^N w_i.
$$

Here:
- $N$ is the number of Monte Carlo samples,
- $w_i$ are the Monte Carlo weights,
- $\bar w$ is the sample mean of the weights,
- $s^2$ is the sample variance of the weights,
- $\mathrm{SE}(\hat I_N)$ is the estimated standard error of the integral.

For sufficiently large $N$, an approximate 95% confidence interval for the
integral is

$$
\hat I_N \pm 1.96\,\mathrm{SE}(\hat I_N).
$$


### Practical tasks

1. Implement the plain Monte Carlo estimator and compute $\hat I_N$ and $\mathrm{SE}(\hat I_N)$ for $N=10^2,10^3,10^4,10^5$ (and optionally $10^6$ if your hardware permits).
2. Plot the estimate $\hat I_N$ with 95% error bars as a function of $N$ (use a log scale for $N$ on the x-axis).
3. Repeat the estimation many times (e.g. 500 independent runs) for fixed $N$ (e.g., $N=10^6$) to build the empirical distribution of $\hat I_N$; plot a histogram and compare its width with the predicted $\mathrm{SE}$.
4. On a log–log plot, show the absolute error $|\hat I_N - I|$  versus $N$. Fit a straight line to verify the slope is approximately $-1/2$.
