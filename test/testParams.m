addpath(genpath('../'))
%%
clearvars
fh=figure;
Nreps=1e1;
cut=90;
%% Test as function of matrix size
subplot(3,1,1)
hold off
clear ae at ae2 at2
Nsamp=1e3;
for M=1:10
Qsqrt=randn(M,M);
Q=Qsqrt*Qsqrt';
Nreps=1e3;
% Estimate:

Qest=nan(M,M,Nreps);
Qest2=nan(M,M,Nreps);
Qtrue=nan(M,M,Nreps);
Qtrue2=nan(M,M,Nreps);
for i=1:Nreps
% Generate data:
X=Qsqrt*randn(size(Q,2),Nsamp);

% Estimate:
Qest(:,:,i)=robCov(X,cut); %My robust estimate
Qest2(:,:,i)=robCov(X,100-2*(100-cut)); %My robust estimate
Qtrue(:,:,i)=X*X'/size(X,2); %Standard, MLE, estimate given a known mean
Qtrue2(:,:,i)=X(:,1:Nsamp*cut/100)*X(:,1:Nsamp*cut/100)'/(Nsamp*cut/100); %Standard, MLE, estimate given a known mean
end

% Visualize
at(M)=mean(sqrt(sum(sum((Q-Qtrue).^2,1),2)),3)/norm(Q,'fro');
at2(M)=mean(sqrt(sum(sum((Q-Qtrue2).^2,1),2)),3)/norm(Q,'fro');
ae(M)=mean(sqrt(sum(sum((Q-Qest).^2,1),2)),3)/norm(Q,'fro');
ae2(M)=mean(sqrt(sum(sum((Q-Qest2).^2,1),2)),3)/norm(Q,'fro');
end
p1=scatter(1:length(at),at,'filled','DisplayName','MLE');
hold on
p3=scatter(1:length(ae),ae2,'filled','DisplayName',['Robust, reject=' num2str(2*(100-cut)) '%']);
p2=scatter(1:length(ae),ae,'filled','DisplayName',['Robust, reject=' num2str(100-cut) '%']);
p1=scatter(1:length(at),at2,'filled','DisplayName',['MLE ' num2str(cut) '% samples']);
title(['Nreps=' num2str(Nreps) ', Nsamples=' num2str(Nsamp) ', no outliers, reject=' num2str(100-cut) '%'])
xlabel('Covariance size (dimension) (M)')
ylabel('Mean |\hat{Q}-Q|_F / |Q|_F')
legend

%% Test as function of data size
subplot(3,1,2)
hold off
clear ae at ae2 ae3
M=6;
for j=1:5
Qsqrt=randn(M,M);
Q=Qsqrt*Qsqrt';
Nreps=1e2;
% Estimate:

Qest=nan(M,M,Nreps);
Qest2=nan(M,M,Nreps);
Qtrue=nan(M,M,Nreps);
for i=1:Nreps
% Generate data:
X=Qsqrt*randn(size(Q,2),10^j);

% Estimate:
Qest(:,:,i)=robCov(X,cut); %My robust estimate
Qest2(:,:,i)=robCov(X,100-2*(100-cut)); %My robust estimate
Qtrue(:,:,i)=X*X'/size(X,2); %Standard, MLE, estimate given a known mean
end

% Visualize
at(j)=mean(sqrt(sum(sum((Q-Qtrue).^2,1),2)),3)/norm(Q,'fro');
ae(j)=mean(sqrt(sum(sum((Q-Qest).^2,1),2)),3)/norm(Q,'fro');
ae2(j)=mean(sqrt(sum(sum((Q-Qest2).^2,1),2)),3)/norm(Q,'fro');
end
p1=scatter(1:length(at),at,'filled','DisplayName','MLE');
hold on
p3=scatter(1:length(ae),ae2,'filled','DisplayName',['Robust, reject=' num2str(2*(100-cut)) '%']);
p2=scatter(1:length(ae),ae,'filled','DisplayName',['Robust, reject=' num2str(100-cut) '%']);
title([num2str(M) ' x ' num2str(M) ' matrix, Nreps=' num2str(Nreps) ', no outliers, reject=' num2str(100-cut) '%'])
xlabel('Log_{10}(sample size)')
ylabel('|\hat{Q}-Q|_F / |Q|_F')
legend
set(gca,'YScale','log')

%% Test with outliers: different iteration number (to show Niter=2 is suff)
subplot(3,1,3)
clear ae at ae2 ae3

M=6;
Nsamp=1e3;
for r=1:11 %Percent outliers
Qsqrt=randn(M,M);
Q=Qsqrt*Qsqrt';
Nreps=1e2;
% Estimate:

Qest=nan(M,M,Nreps);
Qest2=nan(M,M,Nreps);
Qest3=nan(M,M,Nreps);
Qtrue=nan(M,M,Nreps);
for i=1:Nreps
% Generate data:
X=Qsqrt*randn(size(Q,2),Nsamp);
X=X+1e1*repmat((rand(1,Nsamp)>(1-r/100)).*(randn(1,Nsamp)),size(Q,2),1); 
%Add r% outliers
% Estimate:
Qest(:,:,i)=robCov(X,cut,1); %My robust estimate
Qest2(:,:,i)=robCov(X,cut,2); %My robust estimate
Qest3(:,:,i)=robCov(X,cut,10); %My robust estimate
Qtrue(:,:,i)=X*X'/size(X,2); %Standard, MLE, estimate given a known mean
end

% Visualize
at(r)=mean(sqrt(sum(sum((Q-Qtrue).^2,1),2)),3)/norm(Q,'fro');
ae(r)=mean(sqrt(sum(sum((Q-Qest).^2,1),2)),3)/norm(Q,'fro');
ae2(r)=mean(sqrt(sum(sum((Q-Qest2).^2,1),2)),3)/norm(Q,'fro');
ae3(r)=mean(sqrt(sum(sum((Q-Qest3).^2,1),2)),3)/norm(Q,'fro');
end
p1=scatter(1:length(at),at,'filled','DisplayName','MLE');
hold on
p3=scatter(1:length(ae),ae,'filled','DisplayName',['Robust, Niter=1']);
p2=scatter(1:length(ae),ae2,'filled','DisplayName',['Robust, Niter=2']);
p4=scatter(1:length(ae),ae3,'filled','DisplayName',['Robust, Niter=10']);
title([num2str(M) ' x ' num2str(M) ' matrix, Nreps=' num2str(Nreps) ', Nsamples=' num2str(Nsamp) ', reject=' num2str(100-cut) '%, single axis noise'])
xlabel('% outliers')
ylabel('|\hat{Q}-Q|_F / |Q|_F')
legend
set(gca,'YScale','log')

%%
savefig(gcf,'testParams.fig')