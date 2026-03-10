# Create polynomial matrix

Create polynomial matrix used in local polynomial regression (see
details).

## Usage

``` r
polynomial_matrix(l, u = abs(l), d0 = 0, d1 = 3)
```

## Arguments

- l, u:

  lower bound (usually negative) and upper bound (usually positive) of
  the polynomial matrix.

- d0, d1:

  lower and polynomial degree of the polynomial matrix.

## Details

`polynomial_matrix()` computes the following matrix \$\$ \begin{pmatrix}
(l)^{d_0} & (l)^{d_0+1} & \cdots&(l)^{d_1}\\ (l+1)^{d_0} & (l+1)^{d_0+1}
& \cdots&(l+1)^{d_1} \\ \vdots & \vdots & \cdots & \vdots \\ (p)^{d_0} &
(p)^{d_0+1} & \cdots&(p)^{d_1} \end{pmatrix} \$\$

## Examples

``` r
# For example to reproduce DAF filters
daf <- lp_filter(horizon = 6, endpoints = "DAF")
#> Error in .jcall("jdplus/filters/base/r/LocalPolynomialFilters", "Ljdplus/toolkit/base/core/math/linearfilters/ISymmetricFiltering;",     "filters", as.integer(horizon), as.integer(degree), kernel,     endpoints, d, tweight, passband): RcallMethod: cannot determine object class
q <- 0
X <- polynomial_matrix(l = -6, u = q, d0 = 0, d1 = 3)
#> Error in .jcall("jdplus/toolkit/base/core/math/linearfilters/LocalPolynomialFilters",     "Ljdplus/toolkit/base/core/math/matrices/FastMatrix;", "z",     .jnull("jdplus/toolkit/base/core/math/matrices/FastMatrix"),     as.integer(l), as.integer(u), as.integer(d0), as.integer(d1)): java.lang.UnsupportedClassVersionError: jdplus/toolkit/base/core/math/linearfilters/IFiniteFilter has been compiled by a more recent version of the Java Runtime (class file version 65.0), this version of the Java Runtime only recognizes class file versions up to 61.0
K <- diag(sapply(-6:q, function(i) get_kernel(horizon = 6)[i]))
#> Error in .jcall("jdplus/toolkit/base/core/data/analysis/DiscreteKernel",     "Ljava/util/function/IntToDoubleFunction;", tolower(kernel),     h): RcallMethod: cannot determine object class
e_1 <- c(1, 0, 0, 0)
q0 <- K %*% X %*% solve(t(X) %*% K %*% X, e_1)
#> Error: object 'K' not found
q0
#> Error: object 'q0' not found
daf[, "q=0"]
#> Error: object 'daf' not found
```
