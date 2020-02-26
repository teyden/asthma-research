# oddsratio 1.0.3 (June 19 2018)

* update functions to work with ggplot2 v3.0.0

# oddsratio 1.0.2 (December 08 2017)

## Minor
  * Add CITATION file

# oddsratio 1.0.0 (June 12 2017)

## Major
  * rename functions (snake_case)

# oddsratio 0.3.1 (Nov 9 2016)

* update functions to work with ggplot2 v2.2.0
* add data and enable lazy loading in examples

# oddsratio 0.3.0 (Oct 27 2016)

#### New functions
* `plot_smooth.gam()`: Lets you plot smoothing functions of GAM(M)s using `ggplot2`.
* `add.oddsratio.into.plot()`: Add odds ratios into plot of GAM(M) smoothing function.

#### Function updates
* `calc.oddsratio.glm`, `calc.oddsratio.gam`: Add odds ratio confident interval calculation 
* For GLM models CI level can be specified manually.
* Print 'CI' warning if model is of type `glmmPQL`

# oddsratio 0.2.0 (Oct 12 2016)

* Remove param `quietly`
* return data.frame in any case
* update DESCRIPTION

# oddsratio 0.1.0 (Oct 11 2016)

* Initial release attempt to CRAN
