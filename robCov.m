function [Q,outliers]=robCov(w,prc,Niter)
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
if nargin<3 || isempty(Niter)
    Niter=2;
end
if nargin<2 || isempty(prc)
    prc=[];%90; %Default: 90%
    %Change to default finds the 'natural' cut-off. Test
    Niter=5; %This requires more iterations
    %fh=figure; hold on; grid on; axis equal; 
elseif prc<1 %Assuming percentile was given in [0,1] range
    prc=round(100*prc);
end


k=getScale(prc,nD,M);
Q=norm(w,'fro').^2/(M*nD)*eye(nD); %MLE spherical estimate
m=[]; %Presuming zero-mean data.

for i=1:Niter
    [p,y]=z2prctile(w,Q,m); %if w~N(m,Q) this is distributed as t^2 ~ Hotelling's T^2 = nD*(M-1)/(M-nD) F_{nD,M-nD}, see https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution
    if ~isempty(prc) %prc level is given
        yPRC=prctile(y,prc);
    else %Auto-choose prc:
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
        k=getScale(prcAuto,nD,M);
    end
    idx=y<yPRC;
    wRob=w(:,idx);
    Q=(wRob*wRob')/(k*M);
end
outliers=~idx;
%Some debugging:
%prc=[.99,.95,.9];
%yPRC=prctile(y,prc*100);
%y99=yPRC(1);
%y95=yPRC(2);
%y90=yPRC(3);
% figure;
% plot(w(1,:),w(2,:),'o'); hold on;
% plot(w(1,y>y90),w(2,y>y90),'mo');
% plot(w(1,y>y95),w(2,y>y95),'go');
% plot(w(1,y>y99),w(2,y>y99),'ro');
% th=0:.1:2*pi;
% x=sin(th);
% y=cos(th);
%
% col={'r','g','m'};
% for j=1:length(prc)
% f99=finv(prc(j),nD,M-nD); %99th percentile of relevant F distribution
% t99=nD*(M-1)/(M-nD) * f99;  %99th of relevant T^2 distribution
% a99=sqrt((t99)./sum([x;y].*(Q\[x;y])));
% plot(a99.*x,a99.*y,col{j})
% end
% %TODO: deal with auto-correlated (ie not white) noise for better estimates.
end

function k=getScale(prc,nD,M)
%First moment of F_{nD,M-nD} = (M-nD)/(M-nD-2)  (~1 if M-nD>>2)
%(approx) Partial first moment to th prcth-percentile:
if prc==100
    k=1; %This needs hardcoding because finv is undefined for 100th percentile
else
X1=finv(prc/100,nD,M-nD);
dx=X1/1e4;
x=[(dx/2):dx:X1];
k=dx*sum(x.*fpdf(x,nD,M-nD)); % k is a factor to account for the usage of only the 'first' 90% samples and still get an ~unbiased estimate.
end
end
