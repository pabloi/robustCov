# robustCov
A simple, heuristic, fast, robust covariance matrix estimator to deal with data with outliers. Implemented in Matlab.
Comparable to Matlab's built-in robustcov() function. robcov() is both more precise in presence of outliers, and about an order of magnitude faster (see test scripts).


To do:
update approach to use EM to classify samples as either outliers or not.
compare EM performance vs. heurisitc where a fixed number of samples are ignored in covariance estimation.
