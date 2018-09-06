function [p,z]=z2prctile(y,Q,m,iQ)
%Percentile corresponding to observations of multivariate normal data y.
%Computed through squared z-scoring. Equivalent to Mahalanobis distance. 
%If y~N(m,Q) (iid), this is distributed as 
%t^2 ~ Hotelling's T^2 = nD*(M-1)/(M-nD) F_{nD,M-nD}
%see https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution
%INPUTS:
%y: matrix or vector of data
%Q: covariance of data, such that y~N(0,Q);
%m: optional, mean of data, such that y~N(m,Q). If not given, defaults to m=0
%iQ: optional, inverse of Q. If given, it is used instead of Q to save a
%matrix inversion computation.
%OUTPUT:
%p: percentile corresponding to the z^2 score of the data
%z: z^2 score of data
%See also: s2score

if nargin<3
    m=[];
end
if nargin<4
    iQ=[];
end
%First, compute z^2 scores:
z=z2score(y,Q,m,iQ);

[nD,M]=size(y);
p=fcdf(z*(M-nD)/(nD*(M-1)),nD,M-nD);

end