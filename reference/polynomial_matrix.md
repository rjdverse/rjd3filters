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
q <- 0
X <- polynomial_matrix(l = -6, u = q, d0 = 0, d1 = 3)
K <- diag(sapply(-6:q, function(i) get_kernel(horizon = 6)[i]))
e_1 <- c(1, 0, 0, 0)
q0 <- K %*% X %*% solve(t(X) %*% K %*% X, e_1)
q0
#>             [,1]
#> [1,] -0.01723665
#> [2,]  0.02188707
#> [3,]  0.04000228
#> [4,] -0.03414681
#> [5,] -0.09789419
#> [6,]  0.13220425
#> [7,]  0.95518406
daf[, "q=0"]
#> [1] " - 0.0172 B^6 + 0.0219 B^5 + 0.0400 B^4 - 0.0341 B^3 - 0.0979 B^2 + 0.1322 B + 0.9552"
```
