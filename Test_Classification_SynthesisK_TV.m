% test using only one labelling function for multi-class  u={1,2,...,K}
clear
clc
rand('seed',3000);
randn('seed',3000);

load GraphDataSyn_4Circles.mat
[p,M]=size(D); 
N=M/4;
Srate=100;
I=randperm(N);Iset4=I(1:Srate);
I=randperm(N);Iset3=N+I(1:Srate);
I=randperm(N);Iset2=2*N+I(1:Srate);
I=randperm(N);Iset1=3*N+I(1:Srate);

K=4;
Iset=[Iset4 Iset3 Iset2 Iset1]; % reverse order for the label

FD0=zeros(M,K); % 
FD0(Iset4,4)=1;
FD0(Iset3,3)=1;
FD0(Iset2,2)=1;
FD0(Iset1,1)=1;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate the graph
% h=1e4;
% k=5;
% [L,d,lambda_max]=GenerateGraph_fun(data,h,K);

G=sparse(diag(d)-L);
tol=1e-5; % Tolerance for ADMM
maxit=500; % Maximum iterations
% initialization of u--
u00=rand(M,K); 
u00(Iset4,4)=1;
u00(Iset3,3)=1;
u00(Iset2,2)=1;
u00(Iset1,1)=1;
u00=projl1p_1D(u00,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FrameType='Haar'; % FrameType='Linear'; % FrameType='Cubic'; % FrameType='Pseudo-Spline31';
[DFilters RFilters]=ExtractMasks(FrameType);
s=2; % Dilation scale
n=10; % n-1 = Degree of Chebyshev Polynomial Approximation
Lev=1; % Level of transform
lambda=1e-2; % As in lambda||Wu||_1
mu=1e-2; % Parameter from ADMM
M=length(L);
J=log(lambda_max/pi)/log(s)+Lev-1; % Dilation level to start the decomposition
W = @(FD)(GraphWFTG_Decomp(FD,L,DFilters,n,s,J,Lev));
WT = @(d)(GraphWFTG_Recon(d,L,RFilters,n,s,J,Lev));

[u1, energy1,residual1,error1]=SplitBregGraphClassK(FD0,Iset,u00,mu,lambda,d,tol,W,WT,maxit,FD);

figure
subplot(131);plot(log10(residual1)),title('Residual (relative)');
subplot(132);plot((energy1)),title('Energy');
subplot(133);plot((error1)),title('Error');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda=0.05; % As in lambda||u||_TV
[u2, energy2,residual2,error2]=TV_PDHGm_ClassK(FD0,Iset,u00,lambda,tol,G,maxit,1,FD);

figure
subplot(131);plot(log10(residual2)),title('Residual (relative)');
subplot(132);plot((energy2)),title('Energy');
subplot(133);plot((error2)),title('Error');
%%

[uc,FDr1]=max(u1,[],2);
FDr1=FDr1-1;
FDr1(Iset)=labels(Iset);
Error1=sum(abs(FDr1-labels))/(M-length(Iset))*100;

[uc,FDr2]=max(u2,[],2);
FDr2=FDr2-1;
FDr2(Iset)=labels(Iset);
Error2=sum(abs(FDr2-labels))/(M-length(Iset))*100;


%%

figure
subplot(131);scatter(D(1,:),D(2,:),5,FD);title('Ground Truth');axis square;
subplot(132);scatter(D(1,:),D(2,:),5,FDr1);title(['WF error = ',num2str(Error1),'%']);axis square;
subplot(133);scatter(D(1,:),D(2,:),5,FDr2);title(['TV  error = ',num2str(Error2),'%']);axis square;



