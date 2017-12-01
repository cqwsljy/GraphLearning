%load('dataset_local_level1_2_3_4_5_squre_smregion_133.mat')
D=A';
FD=B;
data=D(1:26,:);%the firt scalar,note the transpose of A
K=2;
[P,M]=size(data); %P is features of the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial parameters about WF model and TV model and ADMM
[L,d,lambda_max]=GenerateGraph_fun(data,20000000,20); % GenerateGraph_fun(data,h,k-nearest)
G=sparse(diag(d)-L);
tol=1e-5; % Tolerance for ADMM
maxit=600; % Maximum iterations
FrameType='Haar'; % FrameType='Linear'; % FrameType='Cubic'; % FrameType='Pseudo-Spline31';
[DFilters, RFilters]=ExtractMasks(FrameType);
s=2; % Dilation scale
n=10; % n-1 = Degree of Chebyshev Polynomial Approximation
Lev=1; % Level of transform
lambda=0.01; % As in lambda||Wu||_1
mu=1e-2; % Parameter from ADMM
M=length(L);
J=log(lambda_max/pi)/log(s)+Lev-1; % Dilation level to start the decomposition
W = @(FD)(GraphWFTG_Decomp(FD,L,DFilters,n,s,J,Lev));
WT = @(d)(GraphWFTG_Recon(d,L,RFilters,n,s,J,Lev));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initilize data and labels for WF and TV model
Na=7; % initialize proportion
u00=zeros(M,K);
FD0=u00;
Iset0=[];
error_TV=zeros(M,1);
error_WF=error_TV;
for j=1:M
    index=1:M;
    index_new=setdiff(index,j);
    for i=0:K-1
        index=find(FD(index_new)==i);
        FD0(index,i+1)=1;
        u00(index,i+1)=1;% initialization of u--
    end
    Iset0=index_new;
    Iset=Iset0;clear Iset0
    
    dataTrain=data(:,index_new);
    labelTrain=FD(index_new);
    SVMStruct = fitcsvm(dataTrain',labelTrain,'Prior','uniform','Standardize',1);
    Prediction = predict(SVMStruct,data(:,j)');
    
%     if Prediction==0
%         FD0(j,1)=1;
%         u00(j,1)=1;
%     else
%         FD0(j,2)=1;
%         u00(j,2)=1;
%     end
    
    
%     [u1, energy1,residual1,error1]=SplitBregGraphClassK(FD0,Iset,u00,mu,lambda,d,tol,W,WT,maxit,FD);%WF
%     [uc1,FDr1]=max(u1,[],2);
%     FDr1=FDr1-1;
%     FDr1(Iset)=FD(Iset);
%     error_WF(j,1)=abs(FDr1(j)-FD(j));
%     
%     
%     
%     [u2, energy2,residual2,error2]=TV_PDHGm_ClassK(FD0,Iset,u00,lambda,tol,G,maxit,1,FD);%TV
%     [uc2,FDr2]=max(u2,[],2);
%     FDr2=FDr2-1;
%     FDr2(Iset)=FD(Iset);
%     error_TV(j,1)=abs(FDr2(j)-FD(j));
end

% save('D:\error_WF','error_WF')
% save('D:\error_TV','error_TV')
