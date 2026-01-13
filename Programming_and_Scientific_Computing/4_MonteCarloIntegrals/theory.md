# Monte Carlo Integration — Theory and Statistical Foundations

## 1. Integrals as expectation values

Let $D \subset \mathbb{R}^d$ be a domain and let $f : D \rightarrow \mathbb{R}$ be an integrable function. We are interested in computing the integral

$$
I = \int_D f(x)\, dx
$$

Let $p(x)$ be some probability density function on $D$ such that

$$
p(x) \ge 0, \quad \int_D p(x)\, dx = 1
$$

> Don't worry about what this probability distribution looks like, for now assume we can create whatever distribution we want!

Recall that from any probability density function, we can write the expectation value of any function $f(x)$ as

$$ \mathbb{E}_p[f(x)] = \int_D f(x) p(x) dx$$

The key here is to convert the original integral into an expectation value. First define a new function $g(x) = f(x)/p(x)$. Our original integral then becomes

$$
I = \int_D g(x)p(x)dx = \mathbb{E}_p[g(x)]
$$

> We are assuming here that $p(x) > 0$ whenever $f(x) \neq 0$. 

Hence, the integral is equal to an expectation value!

Ok, you may be asking yourself, why does it matter? This is meaningful because calculating an expectation value, or more precisely approximating it, is an easier task than computing integrals. 

Therefore, the next step is to use cleaver statistics to estimate $\mathbb{E}_p[g(x)]$.

## 2. Monte Carlo estimator

Let $X_1, X_2, \dots, X_N$ be independent and identically distributed random variables drawn from the probability density $p(x)$.

The expectation value above is approximated by the empirical mean

$$
\mathbb{E}_p\!\left[\frac{f(X)}{p(X)}\right]
\approx
\frac{1}{N} \sum_{i=1}^N \frac{f(X_i)}{p(X_i)}
$$

which defines the Monte Carlo estimator

$$
\hat{I}_N
=
\frac{1}{N} \sum_{i=1}^N \frac{f(X_i)}{p(X_i)}
$$


## 4. Statistical error and the Central Limit Theorem

Assume that the variance

$$
\sigma^2
=
\mathrm{Var}_p\!\left[\frac{f(X)}{p(X)}\right]
$$

is finite.

Then, by the central limit theorem,

$$
\sqrt{N}(\hat{I}_N - I)
\;\xrightarrow{d}\;
\mathcal{N}(0, \sigma^2)
$$

which implies that for large $N$

$$
\hat{I}_N \approx \mathcal{N}
\left(
I,
\frac{\sigma^2}{N}
\right)
$$

This result explains the universal $1/\sqrt{N}$ convergence rate of Monte Carlo integration.

---

## 5. Error estimation from samples

The variance $\sigma^2$ is typically unknown but can be estimated from the samples.

Define

$$
w_i = \frac{f(X_i)}{p(X_i)}
$$

and the sample variance

$$
s^2 = \frac{1}{N-1}
\sum_{i=1}^N
\left(
w_i - \frac{1}{N}\sum_{j=1}^N w_j
\right)^2
$$

An estimate of the standard error of the Monte Carlo integral is

$$
\Delta I \approx \frac{s}{\sqrt{N}}
$$

---

## 6. Choice of sampling distribution

The choice of the sampling distribution $p(x)$ directly affects the variance of the estimator.

The optimal distribution (in the sense of minimal variance) is

$$
p^\star(x) \propto |f(x)|
$$

for which the variance vanishes formally. In practice, this distribution is rarely accessible, but it motivates the design of **importance sampling** strategies where $p(x)$ approximates the shape of $|f(x)|$.

Poor choices of $p(x)$, especially those that are small where $|f(x)|$ is large, can lead to very large variances or even divergence.


## 7. Special case: uniform sampling

If $p(x)$ is uniform over a bounded domain $D$,

$$
p(x) = \frac{1}{|D|}
$$

then the general estimator reduces to

$$
\hat{I}_N
=
\frac{|D|}{N} \sum_{i=1}^N f(X_i)
$$

recovering the plain Monte Carlo estimator introduced in elementary examples such as the estimation of $\pi$.


## 8. High-dimensional integrals

The convergence rate of Monte Carlo integration does not depend on the dimension $d$ of the integration domain. While the variance may increase with dimension due to the behavior of $f(x)$, the asymptotic error still scales as $1/\sqrt{N}$.

This property makes Monte Carlo integration particularly suitable for high-dimensional problems, where deterministic quadrature methods suffer from exponential complexity.
