# SBC and various test quantities

Simulations for the paper "Simulation-Based Calibration Checking for Bayesian Computation: The Choice of Test Quantities Shapes Sensitivity". Written in
R markdown, using the [SBC](https://hyunjimoon.github.io/SBC/) R package.
Code is written in R markdown (`.Rmd`), rendered HTML can be found at
https://martinmodrak.github.io/sbc_test_quantities_paper.

Contents:

- `mvn.Rmd` - the examples for the main body of the paper (Section 4 - Numerical Experiments)
- `bernoulli.Rmd` - Examples in the Appendix B

The code uses `renv` to recreate the same environment as we used to run the examples.
Run `renv::restore()` after loading the project.

The rendered site can be rebuilt with `rmarkdown::render_site()`.
