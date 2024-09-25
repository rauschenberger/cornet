
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/rauschenberger/cornet?svg=true)](https://ci.appveyor.com/project/rauschenberger/cornet)
[![R-CMD-check](https://github.com/rauschenberger/cornet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rauschenberger/cornet/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/rauschenberger/cornet/graph/badge.svg)](https://app.codecov.io/gh/rauschenberger/cornet)

# Penalised regression for dichotomised outcomes 

Implements lasso and ridge regression for dichotomised outcomes (i.e., numerical outcomes that were transformed to binary outcomes).

## Installation

Install the current release from
[CRAN](https://CRAN.R-project.org/package=cornet):

``` r
install.packages("cornet")
```

or the latest development version from
[GitHub](https://github.com/rauschenberger/cornet):

``` r
#install.packages("remotes")
remotes::install_github("rauschenberger/cornet")
```

## Reference

Armin Rauschenberger
[![AR](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0001-6498-4801)
and Enrico Glaab
[![EG](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-3977-7469)
(2024). "Predicting dichotomised outcomes from high-dimensional data in biomedicine".
*Journal of Applied Statistics* 51(9):1756-1771.
[doi: 10.1080/02664763.2023.2233057](https://doi.org/10.1080/02664763.2023.2233057).

[![CRAN version](https://www.r-pkg.org/badges/version/cornet)](https://CRAN.R-project.org/package=cornet)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/cornet)](https://CRAN.R-project.org/package=cornet)
[![Total CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/cornet)](https://CRAN.R-project.org/package=cornet)

## Disclaimer

The R package `cornet` implements elastic net regression for dichotomised outcomes ([Rauschenberger & Glaab, 2024](https://doi.org/10.1080/02664763.2023.2233057)).

Copyright &copy; 2018 Armin Rauschenberger, University of Luxembourg, Luxembourg Centre for Systems Biomedicine (LCSB), Biomedical Data Science (BDS)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
