% test using only one labelling function for multi-class  u={1,2,...,K}
clear
clc

rand('seed',3000);
randn('seed',3000);

load minst_5e3.mat;
[p,M]=size(data);
Srate=round(M*0.05);
I=randperm(M);Iset=I(1:Srate);
clear I;
FD0=zeros(M,K); % 
FD0(Iset,labels(Iset)+1)=1;

FD=labels;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate the graph
h=1e4;
k=5;
[L,d,lambda_max]=GenerateGraph_fun(data,h,K); % Laplacian matrix
G=sparse(diag(d)-L); % Weighted adjacent matrix


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=10;
tol=1e-5; % Tolerance for ADMM
maxit=500; % Maximum iterations

% intialization of u00
[IDX,c0] = kmeans(data,K);
u00=zeros(M,K);
 for k=1:K
       idx=(IDX==k);
   u00(idx,k)=1;
 end
 u00(Iset,labels(Iset)+1)=1;
 %%
%u00=projl1p_1D(rand(M,K),1);

%% parameter set for Wavelet frame based
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

%
[u1, energy1,residual1,error1]=SplitBregGraphClassK(FD0,Iset,u00,mu,lambda,d,tol,W,WT,maxit,FD);

%%
lambda=1e-2; % As in lambda||u||_TV
[u2, energy2,residual2,error2]=TV_PDHGm_ClassK(FD0,Iset,u00,lambda,tol,G,maxit,1,FD);

%% Compute the final classification and rate
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
subplot(221);plot(log10(energy1));title('WF energy');
subplot(222);plot(log10(energy2));title('TV energy');
subplot(223);plot((error1));title(['Wavelet frame: error = ',num2str(Error1),'%']);axis square;
subplot(224);plot((error2));title(['TV: error = ',num2str(Error2),'%']);axis square;