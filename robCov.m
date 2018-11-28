function [Q,outliers,prc]=robCov(w,prc,Niter)
%robCov is a robust covariance matrix estimation of data. Useful if data may contain outliers.
%It uses only an 'inner' percentage of the data (i.e. data lying inside a certain
%ellipsoid) to estimate the covariance matrix. If data comes from a
%multinormal distribution, this estimation is unbiased. There is no
%guarantee for data coming from other distributions (e.g. heavier tailed ones).
%Data is assumed to be 0-mean.
%The procedure is as follows:
%1) Estimate the covariances with the classical MLE estimator.
%2) Use the classical estimator to find the X% of inner-most samples (i.e.
%samples within the 90-th percentile ellipsoid). The other (100-X)% is
%presumed to contain all the outliers.
%3) Estimate the covariance of the reduced samples.
%4) Expand estimate to account for the fact that only the samples closes to
%the origin where used, under the assumption of multivariate normal data.
%INPUT:
%w: MxN matrix, consisting of N-samples of M-dimensional data.
%prc: (optional, default =90%) Percentage of data used for estimation.
%OUTPUT:
%Q: MxM covariance estimate
%See also: robustCov, estimateParams

%To do:
%1) is it worth it to iterate the procedure to converge on the outlier samples?
%2) Instead of hard-coding the % of samples to discard, do auto-setting:
%compute the z^2 scores percentiles, and set the cut-off whenever the
%empirical percentiles are outside some window of the expected percentiles
%(e.g. if the sample datapoint at the 95% percentile is given a percentile score of more than 96% or less than 94%, cut off there)

[nD,M]=size(w);
%Init:
Q=norm(w,'fro').^2/(M*nD)*eye(nD); %MLE spherical estimate
m=[]; %Presuming zero-mean data.

%Parse inputs:
if nargin<3 || isempty(Niter)
    Niter=2;
end
if nargin<2 || isempty(prc) %AUTO-prc mode
    prc=[];%90; %Default: 90%
    %Change to default finds the 'natural' cut-off. Test
    Niter=5; %This requires more iterations
    %fh=figure; hold on; grid on; axis equal; 
else
    if prc<1 %Assuming percentile was given in [0,1] range
        prc=round(100*prc);
    end
    M2=round(M*prc/100); 
    k=getScale(prc,nD,M2);
end

%Q=norm(w,'fro').^2/(M*nD)*eye(nD); %MLE spherical estimate
w2=w*w'; 
Q=(w2)/M; %Standard estimate, to init 
m=[]; %Presuming zero-mean data.


for i=1:Niter
    %Q
    if ~isempty(prc) %prc level is given
        y=z2score(w,Q,m); %if w~N(m,Q) this is distributed as t^2 ~ Hotelling's T^2 = nD*(M-1)/(M-nD) F_{nD,M-nD}, see https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution
        yPRC=prctile(y,prc);
    else %Auto-choose prc:
        [p,y]=z2prctile(w,Q,m,[],M); %Computing z-scores and associated percentiles given the underlying distribution, should be done with M2 if those are the samples used for estimating Q
        pEmp=[1:length(p)]/length(p);
        aux=(sort(p)-pEmp); %Divergence between empirical (ranking) percentiles, and expected percentiles given the latest covariance estimate
        prcOut1=max(aux); %Peak
        prcOut2=1-find(aux==max(aux),1,'last')/M; %Where the peak ocurrs
        prcOut=2*min(prcOut1,prcOut2); %Taking minimum estimate, and adding safety margin
        if prcOut<0 || isempty(prcOut)
            prcOut=.001; %Defaults to .1% outliers if nothing was found
        end
        prcAuto=max(100*(1-prcOut),60); %No more than 40% outliers
        %plot([1:M]/M,aux); plot(1-prcOut1,0,'kx'); plot(1-prcOut2,0,'ko'); plot(prcAuto/100,0,'rd'); drawnow;
        yPRC=prctile(y,prcAuto);
        M2=round(M*prcAuto/100); 
        k=getScale(prcAuto,nD,M2);
    end
    idx=y<yPRC; %Non-rejection criteria
    wRob=w(:,idx);
    Q=(wRob*wRob')/(k*M); %Should I use ACTUAL DoF instead of expected DoF?
end
outliers=~idx;
if isempty(prc)
    prc=prcAuto;
end
end

function k=getScale(prc,nD,M,fc)
%First moment of F_{nD,M-nD} = (M-nD)/(M-nD-2)  (~1 if M-nD>>2)
%(approx) Partial first moment to th prcth-percentile:
if prc>99.9 %prc is essentially 100%
    k=fstat(nD,M-nD); %This needs hardcoding to avoid sampling fpdf in too many points
else
    if nargin<4 %Compute everything
        X1=finv(prc/100,nD,M-nD);
        dx=1e-3;
        x=[(dx/2):dx:X1];
        fp=fpdf(x,nD,M-nD);
        k=dx*sum(x.*fp); 
    else %fc was precomputed
        X1=finv(prc/100,nD,M-nD);
        Xidx=round(1e3*X1); %Round to 1e-3
        k=fc(Xidx);
    end
end
end
