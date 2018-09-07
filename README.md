# robustCov
A simple, heuristic, fast, robust covariance matrix estimator to deal with data with outliers. Implemented in Matlab.
Comparable to Matlab's built-in robustcov() function. robcov() is both more precise in presence of outliers, and about an order of magnitude faster (see test scripts).

%TO DO:
- Correct expansion factor by estimating true number of outliers after covariance has been estimated, as % data outside the 99% CI (minus 1%).
- Debug iterative approach that leads to performance downgrading for Niter>2
- Improve automatic percentile selection: optimize for Niter and safety margin parameters.
