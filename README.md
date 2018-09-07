# robustCov
A simple, fast, robust covariance matrix estimator to deal with data with outliers. Implemented in Matlab.

%TO DO:
- Correct expansion factor by estimating true number of outliers after covariance has been estimated, as % data outside the 99% CI (minus 1%).
- Debug iterative approach that leads to performance downgrading for Niter>2
- Improve automatic percentile selection: optimize for Niter and safety margin parameters.
