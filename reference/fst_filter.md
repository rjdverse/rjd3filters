# Estimation of a filter using the Fidelity-Smoothness-Timeliness criteria

Estimation of a filter using the Fidelity-Smoothness-Timeliness criteria

## Usage

``` r
fst_filter(
  lags = 6,
  leads = 0,
  pdegree = 2,
  smoothness.weight = 1,
  smoothness.degree = 3,
  timeliness.weight = 0,
  timeliness.passband = pi/6,
  timeliness.antiphase = TRUE
)
```

## Arguments

- lags:

  Lags of the filter (should be positive).

- leads:

  Leads of the filter (should be positive or 0).

- pdegree:

  Local polynomials preservation: max degree.

- smoothness.weight:

  Weight for the smoothness criterion (in \\\[0, 1\]\\).

- smoothness.degree:

  Degree of the smoothness criterion (3 for Henderson).

- timeliness.weight:

  Weight for the Timeliness criterion (in \\\[0, 1\[\\).
  `sweight+tweight` should be in \\\[0,1\]\\.

- timeliness.passband:

  Passband for the timeliness criterion (in radians). The phase effect
  is computed in \\\[0, passband\]\\.

- timeliness.antiphase:

  boolean indicating if the timeliness should be computed analytically
  (`TRUE`) or numerically (`FALSE`).

## Details

Moving average computed by a minimisation of a weighted mean of three
criteria under polynomials constraints. Let \\\boldsymbol
\theta=(\theta\_{-p},\dots,\theta\_{f})'\\ be a moving average where
\\p\\ and \\f\\ are two integers defined by the parameter `lags` and
`leads`. The three criteria are:

- *Fidelity*, \\F_g\\: it's the variance reduction ratio. \$\$
  F_g(\boldsymbol \theta) = \sum\_{k=-p}^{+f}\theta\_{k}^{2} \$\$

- *Smoothness*, \\S_g\\: it measures the flexibility of the coefficient
  curve of a filter and the smoothness of the trend. \$\$
  S_g(\boldsymbol \theta) = \sum\_{j}(\nabla^{q}\theta\_{j})^{2} \$\$
  The integer \\q\\ is defined by parameter `smoothness.degree`. By
  default, the Henderson criteria is used (`smoothness.degree = 3`).

- *Timeliness*, \\T_g\\ : \$\$
  T_g(\boldsymbol\theta)=\int\_{0}^{\omega\_{2}}f(\rho\_{\boldsymbol\theta}(\omega),\varphi\_{\boldsymbol\theta}(\omega))d\omega
  \$\$ with \\\rho\_{\boldsymbol\theta}\\ and
  \\\varphi\_{\boldsymbol\theta}\\ the gain and phase shift functions of
  \\\boldsymbol \theta\\, and \\f\\ a penalty function defined as
  \\f\colon(\rho,\varphi)\mapsto\rho^2\sin(\varphi)^2\\ to have an
  analytically solvable criterium. \\\omega\_{2}\\ is defined by the
  parameter `timeliness.passband` and is it by default equal to
  \\2\pi/12\\: for monthly time series, we focus on the timeliness
  associated to cycles of 12 months or more.

The moving average is then computed solving the problem: \$\$
\begin{cases} \underset{\theta}{\min} & J(\theta)= (1-\beta-\gamma)
F_g(\theta)+\beta S_g(\theta)+\gamma T_g(\theta)\\ s.t. & C\theta=a
\end{cases} \$\$ Where \\C\theta=a\\ represents linear constraints to
have a moving average that preserve polynomials of degree \\q\\
(`pdegree`): \$\$ C=\begin{pmatrix} 1 & \cdots&1\\ -h & \cdots&h \\
\vdots & \cdots & \vdots \\ (-h)^d & \cdots&h^d \end{pmatrix},\quad
a=\begin{pmatrix} 1 \\0 \\ \vdots\\0 \end{pmatrix} \$\$

## References

Grun-Rehomme, Michel, Fabien Guggemos, and Dominique Ladiray (2018).
“Asymmetric Moving Averages Minimizing Phase Shift”. In: Handbook on
Seasonal Adjustment,
<https://ec.europa.eu/eurostat/web/products-manuals-and-guidelines/-/ks-gq-18-001>.

## Examples

``` r
filter <- fst_filter(lags = 6, leads = 0)
filter
#> [1] "0.1678 B^6 - 0.3147 B^4 - 0.2797 B^3 + 0.2098 B^2 + 0.6713 B + 0.5455"
```
