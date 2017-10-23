%TV and WF comparasiom code
% load('dataset_local_level1_2_3_4_5_squre_smregion_133.mat')
load('E:\MultiscaleAD_Classification\ADNI_screen1.5T_all_level6_mwp1_hippo.mat')
% D=A';
% FD=B;
%D=D(1:52,:);%the firt scalar of wavelete,note the transpose of A
FD=FD-1;
K=4;
[P,M]=size(D);



[L,d,lambda_max]=GenerateGraph_fun(D,200000,20); % GenerateGraph_fun(data,h,k-nearest)

G=sparse(diag(d)-L);
tol=1e-5; % Tolerance for ADMM
maxit=600; % Maximum iterations

FrameType='Haar'; % FrameType='Linear'; % FrameType='Cubic'; % FrameType='Pseudo-Spline31';
[DFilters RFilters]=ExtractMasks(FrameType);
s=2; % Dilation scale
n=10; % n-1 = Degree of Chebyshev Polynomial Approximation
Lev=1; % Level of transform
lambda=0.01; % As in lambda||Wu||_1
mu=1e-2; % Parameter from ADMM
M=length(L);
J=log(lambda_max/pi)/log(s)+Lev-1; % Dilation level to start the decomposition
W = @(FD)(GraphWFTG_Decomp(FD,L,DFilters,n,s,J,Lev));
WT = @(d)(GraphWFTG_Recon(d,L,RFilters,n,s,J,Lev));
iters=1;
err=zeros(iters,2);
Na=10;%each class has about 71 samples,Na/71 is proportion
for j=1:iters  
    u00=zeros(M,K);
    FD0=u00;
    Iset0=[];
    for i=0:K-1
        index=find(FD==i);
        index=index(randperm(length(index)));
        index=index(1:Na);
        FD0(index,i+1)=1;
        Iset0=[Iset0;index];
        u00(index,i+1)=1;% initialization of u--
    end
    Iset=Iset0;clear Iset0
    %%%SplitBregGraphClassK
    [u1, energy1,residual1,error1]=SplitBregGraphClassK(FD0,Iset,u00,mu,lambda,d,tol,W,WT,maxit,FD);
    %{
    figure
    subplot(131);plot(log10(residual1)),title('WF Residual (relative)');
    subplot(132);plot(log10(energy1)),title('WF Energy');
    subplot(133);plot((error1)),title('WF Error');
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    lambda=0.05; % As in lambda||u||_TV
    %%%TV_PDHGm_ClassK
    [u2, energy2,residual2,error2]=TV_PDHGm_ClassK(FD0,Iset,u00,lambda,tol,G,maxit,1,FD);
    %{
    figure
    subplot(131);plot(log10(residual2)),title('TV Residual (relative)');
    subplot(132);plot(log10(energy2)),title('TV Energy');
    subplot(133);plot((error2)),title('TV Error');
    %%
    %}
    [uc,FDr1]=max(u1,[],2);
    FDr1=FDr1-1;
    FDr1(Iset)=FD(Iset);
    Error1=sum(abs(FDr1-FD))/(M-length(Iset))*100;
    
    [uc,FDr2]=max(u2,[],2);
    FDr2=FDr2-1;
    FDr2(Iset)=FD(Iset);
    Error2=sum(abs(FDr2-FD))/(M-length(Iset))*100;
    %
    %{
    figure
    subplot(131);scatter(D(1,:),D(2,:),5,FD);title('Ground Truth');axis square;
    subplot(132);scatter(D(1,:),D(2,:),5,FDr1);title(['WF error = ',num2str(Error1),'%']);axis square;
    subplot(133);scatter(D(1,:),D(2,:),5,FDr2);title(['TV  error = ',num2str(Error2),'%']);axis square;
    %}
    err(j,1)=Error1;
    err(j,2)=Error2;
end
[~,index]=sort(FD);
FD=FD(index);
D=D(:,index);
FD=FD(1:412);
D=D(:,1:412);
[L]=GraphLaplacian(D,2e10,20);%standard spectral
[vector,value]=eig(L);
[~,idx]=sort(abs(diag(value)));
y=vector(:,idx(1:K));
label=kmeans(y,K);
label=label-1;
error5=abs(label-FD);%error5 is a vector
% if sum(error5)/size(data,2)>0.5
%     error5=abs(error5-1);
% end
% Error5=100*sum(error5)/133;

