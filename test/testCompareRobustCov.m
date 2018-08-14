%% Compare performance and runtime to robustCov()
fh=figure;
clear ae at ae2 ae3 rcTime rc2Time time
M=6;
Nsamp=1e3;
for r=1:5 %Percent outliers
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
X=X+1e1*repmat((rand(1,Nsamp)>(1-2*r/100)).*(randn(1,Nsamp)),size(Q,2),1); 
%Add r% outliers
% Estimate:
tic
Qest(:,:,i)=robCov(X,cut); %My robust estimate
rcTime(i)=toc;
tic
Qest2(:,:,i)=robustcov(X'); %MAtlab built-in
rc2Time(i)=toc;
tic
Qtrue(:,:,i)=X*X'/size(X,2); %Standard, MLE, estimate given a known mean
time(i)=toc;
end

% Visualize
at(r)=sqrt(mean(sum(sum((Q-Qtrue).^2,1),2),3))/norm(Q,'fro');
ae(r)=sqrt(mean(sum(sum((Q-Qest).^2,1),2),3))/norm(Q,'fro');
ae2(r)=sqrt(mean(sum(sum((Q-Qest2).^2,1),2),3))/norm(Q,'fro');
rt(r)=mean(time);
re(r)=mean(rcTime);
re2(r)=mean(rc2Time);
end
subplot(2,1,1)
p1=scatter([1:length(at)]*2,at,'filled','DisplayName','MLE');
hold on
p3=scatter([1:length(at)]*2,ae,'filled','DisplayName',['robCov()']);
p2=scatter([1:length(at)]*2,ae2,'filled','DisplayName',['robustcov()']);
title(['Performance comparison to robustcov() [same conds as above]'])
xlabel('% outliers')
ylabel('|\hat{Q}-Q|_F / |Q|_F')
legend
set(gca,'YScale','log')
subplot(2,1,2)
p1=scatter([1:length(at)]*2,rt,'filled','DisplayName','MLE');
hold on
p3=scatter([1:length(at)]*2,re,'filled','DisplayName',['RobCov()']);
p2=scatter([1:length(at)]*2,re2,'filled','DisplayName',['RobustCov()']);
xlabel('% outliers')
ylabel('Avg. run time (s)')
legend
set(gca,'YScale','log')
title('Run time comparison to robustcov()')