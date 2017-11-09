% test using only one labelling function for multi-class  u={1,2,...,K}
clear
clc

% load GraphDataSyn_4Circles.mat
load GraphDataSyn_4CirclesZ2_5nn.mat
[p,M] = size(D); %
classK = length(unique(FD));
proInitial = 0.0357;
Srate = round(proInitial * M /classK);

FD0 = zeros(M,classK);
u00 = zeros(M,classK);
Iset = [];
for k = 1:classK
	index = find(FD == k - 1);
	I = randperm(length(index));
	index = index(I);
	Isetk = ['Iset',num2str(k)];
	eval('Isetk=index(1:Srate);');
	eval('FD0(Isetk,k) = 1;');
	eval('u00(Isetk,k) = 1;');
    eval('Iset = [Iset;Isetk];');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate the graph
h = 1e4;
Knears = 5;
% [L,d,lambda_max]=GenerateGraph_fun(D,h,Knears,'ZM2'); 

sindex = 1:length(d);
G = sparse(sindex,sindex,d)-L; %é‚»æŽ¥çŸ©é˜µ
tol = 1e-5; % Tolerance for ADMM
maxit = 500; % Maximum iterations
clear sindex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FrameType = 'Haar'; % FrameType='Linear'; % FrameType='Cubic'; % FrameType='Pseudo-Spline31';
[DFilters, RFilters] = ExtractMasks(FrameType);
s = 2; % Dilation scale
n = 10; % n-1 = Degree of Chebyshev Polynomial Approximation
Lev = 1; % Level of transform
lambda = 0.01; % As in lambda||Wu||_1
mu = 1e-1; % Parameter from ADMM
M = length(L);
J = log(lambda_max/pi)/log(s)+Lev-1; % Dilation level to start the decomposition
W = @(FD)(GraphWFTG_Decomp(FD,L,DFilters,n,s,J,Lev));
WT = @(d)(GraphWFTG_Recon(d,L,RFilters,n,s,J,Lev));
disp('WF Model...')
[u1, energy1,residual1,error1] = SplitBregGraphClassK(FD0,Iset,u00,mu,lambda,d,tol,W,WT,maxit,FD);
%[u1, energy1,residual1,error1] = WF_PDHGm_ClassK(FD0,Iset,u00,lambda,d,tol,W,WT,maxit,1,FD);
figure
subplot(131);plot(log10(residual1),'LineWidth',2),title('WF Residual (relative)');
xlabel('Iteration')
ylabel('Log Of Residual')
subplot(132);plot(log10(energy1),'LineWidth',2),title('WF Energy');
xlabel('Iteration')
ylabel('Log Of Energy')
subplot(133);plot((error1),'LineWidth',2),title('WF Error');
xlabel('Iteration')
ylabel('Error')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = 0.05; % As in lambda||u||_TV
disp('TV Model...')
[u2, energy2,residual2,error2] = TV_PDHGm_ClassK(FD0,Iset,u00,lambda,tol,G,maxit,1,FD);


figure
subplot(131);plot(log10(residual2),'LineWidth',2),title('TV Residual (relative)');
xlabel('Iteration')
ylabel('Log Of Residual')
subplot(132);plot(log10(energy2),'LineWidth',2),title('TV Energy');
xlabel('Iteration')
ylabel('Log Of Energy')
subplot(133);plot((error2),'LineWidth',2),title('TV Error');
xlabel('Iteration')
ylabel('Error')
%%

[~,FDr1] = max(u1,[],2);
FDr1 = FDr1-1;
FDr1(Iset) = FD(Iset);
Error1 = sum(abs(FDr1-FD))/(M-length(Iset))*100;

[uc,FDr2] = max(u2,[],2);
FDr2 = FDr2-1;
FDr2(Iset) = FD(Iset);
Error2 = sum(abs(FDr2-FD))/(M-length(Iset))*100;
%
figure
subplot(131);scatter(D(1,:),D(2,:),5,FD);title('Ground Truth');axis square;
subplot(132);scatter(D(1,:),D(2,:),5,FDr1);title(['WF error = ',num2str(Error1),'%']);axis square;
subplot(133);scatter(D(1,:),D(2,:),5,FDr2);title(['TV  error = ',num2str(Error2),'%']);axis square;



