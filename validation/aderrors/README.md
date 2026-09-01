# Gamma Method Validation (Python vs. Julia implementation by Alberto Ramos) 
The Gamma-method implemented in Python is compared to the Julia implementation of the Gamma-method by Alberto Ramos. Alberto's code and the corresponding documentation can be found here:
https://ific.uv.es/~alramos/software/aderrors/

## Test dataset

- Observable: plaquette
- beta: 6
- lattice: 8x8x8x8
- trajectories: 1000
- trajectory length: 1
- molecular dynamics steps: 10
- thermalization cut: 100
- measurements: 910 (10 warmup trajectories)

## Input

`input/plaquette_beta6.csv`

The input data contains the same thermalized time series
used by the Python analysis.

## Quantities compared

- mean
- integrated autocorrelation time
- optimal window
- statistical error