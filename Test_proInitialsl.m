%% test using only one labelling function for multi-class  u={1,2,...,K}
% clear
% clc
% load HippoMCIcvsMCInc133.mat
% load HippoWaveleteAD_NC_MCI_811.mat
% load HippoWaveleteAD_NC_416.mat

% load mnistAll.mat
% load mnistAllZM1.mat
% load('mnist49Z2.mat')
% load mnist49Z2.mat

% load GraphDataSyn_4CirclesZ2_5nn.mat
% load('E:\data\GraphDataSyn_3Circles.mat')

% generate the graph
% h=1e4;
% Knears = 10;
% [L,d,lambda_max]=GenerateGraph_fun(D,h,Knears,'ZM2'); 

[p,M] = size(D); %
classK = length(unique(FD));
proInitials = [0.01,0.03,0.1,0.15,0.2];
errors = zeros(length(proInitials),3);
pro = 1;
for proInitial = proInitials
	    
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
	    eval('Iset = [Iset,Isetk];');
	end
	tol = 1e-5; % Tolerance
	maxit = 600; % Maximum iterations
	%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	sindex = 1:length(d);
	G = sparse(sindex,sindex,d)-L; %
	% clear sindex

	FrameType = 'Linear'; %Haar, FrameType='Linear'; % FrameType='Cubic'; % FrameType='Pseudo-Spline31';
	[DFilters, RFilters] = ExtractMasks(FrameType);
	s = 2; % Dilation scale
	n = 10; % n-1 = Degree of Chebyshev Polynomial Approximation
	Lev = 1; % Level of transform
	lambda = 0.01; % As in lambda||Wu||_1
	mu = 1e-2; % Parameter from ADMM
	M = length(L);
	J = log(lambda_max/pi)/log(s)+Lev-1; % Dilation level to start the decomposition
	W = @(FD)(GraphWFTG_Decomp(FD,L,DFilters,n,s,J,Lev));
	WT = @(d)(GraphWFTG_Recon(d,L,RFilters,n,s,J,Lev));

	disp('WF Model By PDHG...')
	[u1,energy1,residual1,error1] = WF_PDHGm_ClassK(FD0,Iset,u00,lambda,d,mu,tol,W,WT,maxit,1,FD);

	disp('WF Model By ADMM...')
	[u2, energy2,residual2,error2] = SplitBregGraphClassK(FD0,Iset,u00,lambda,d,mu,tol,W,WT,maxit,FD);

	disp('WF Model By Adaptive PDHG...')
	[u5, energy5,residual5,error5] = WFPDHGmApativeClassK(FD0,Iset,u00,lambda,d,mu,tol,W,WT,maxit,FD);
	errors(pro,:) = [100 - min(error1),100 - min(error2),100 - min(error5)];
	pro = pro + 1;
end
%{
figure
subplot(131);plot(log10(residual1)),title('WF Residual (relative)');
subplot(132);plot(log10(energy1)),title('WF Energy');
subplot(133);plot((error1)),title('WF Error');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



lambda = 0.05; % As in lambda||u||_TV
disp('TV Model By PDHG...')
[u3, energy3,residual3,error3] = TV_PDHGm_ClassK(FD0,Iset,u00,lambda,tol,G,maxit,1,FD);

disp('TV Model By ADMM...')
[u4, energy4,residual4,error4] = TV_SplitBregClassK(FD0,Iset,u00,tol,G,maxit,mu,FD);


figure
subplot(131);plot(log10(residual3)),title('TV Residual (relative)');
subplot(132);plot(log10(energy3)),title('TV Energy');
subplot(133);plot((error3)),title('TV Error');
%%


FDr1(Iset) = FD(Iset);
Error1 = sum(abs(FDr1-FD))/(M-length(Iset))*100;

[uc,FDr2] = max(u2,[],2);
FDr2 = FDr2-1;
FDr2(Iset) = FD(Iset);
Error2 = sum(abs(FDr2-FD))/(M-length(Iset))*100;
%

% figure
% subplot(131);scatter(D(1,:),D(2,:),5,FD);title('Ground Truth');axis square;
% subplot(132);scatter(D(1,:),D(2,:),5,FDr1);title(['WF error = ',num2str(Error1),'%']);axis square;
% subplot(133);scatter(D(1,:),D(2,:),5,FDr2);title(['TV  error = ',num2str(Error2),'%']);axis square;
% img1 = [];
% img2 = [];
% for k = 0:9
%     c = FD == k;
%     index = find(c);
%     if k <=4
%         t = D(:,index(100));
%         t = reshape(t,28,28)';
%         img1 = [img1,t];
%     else
%         t = D(:,index(100));
%         t = reshape(t,28,28)';
%         img2 = [img2,t];
%     end
% end
% img = [img1;img2];

% 小波系数画图
% subplot(131);plot(FD,'LineWidth',1.5)
% xlabel('x')
% ylabel('u')
% title('原始标签')
% subplot(132);plot(y{2,1})
% xlabel('x')
% ylabel('Wu高频部分')
% title('小波系数')
% subplot(133);plot(y{2,1})
% xlabel('x')
% ylabel('Wu高频部分')
% title('部分小波系数')
%}
%{

%% plot PDHG and ADMM
subplot(121);
plot(energy1,'LineWidth',1.5);hold on
plot(energy2,'LineWidth',1.5)
legend('PDHG','ADMM')
xlabel('iteration')
ylabel('energy')

subplot(122);
plot(error1,'LineWidth',1.5);hold on
plot(error2,'LineWidth',1.5)
xlabel('iteration')
ylabel('error')
legend('PDHG','ADMM')

%% plot Adaptive PDHG and PDHG and ADMM
subplot(121);
plot(energy1,'LineWidth',1.5);hold on
plot(energy2,'LineWidth',1.5);
plot(energy5,'LineWidth',1.5)
legend('Adaptive PDHG','ADMM','PDHG')
xlabel('iteration')
ylabel('energy')

subplot(122);
plot(error1,'LineWidth',1.5);hold on
plot(error2,'LineWidth',1.5)
plot(error5,'LineWidth',1.5)
xlabel('iteration')
ylabel('error')
legend('Adaptive PDHG','ADMM','PDHG')

%% plot haar and linear
subplot(231);
plot(energy11,'LineWidth',1.5);hold on
plot(energy1,'LineWidth',1.5)
legend('Haar','Linear')
xlabel('iteration')
ylabel('energy')
title('WF PDHG')

subplot(232);
plot(energy21,'LineWidth',1.5);hold on
plot(energy2,'LineWidth',1.5)
legend('Haar','Linear')
xlabel('iteration')
ylabel('energy')
title('WF ADMM')

subplot(233);
plot(energy51,'LineWidth',1.5);hold on
plot(energy5,'LineWidth',1.5)
legend('Haar','Linear')
xlabel('iteration')
ylabel('energy')
title('WF Adaptive PDHG')


subplot(234);
plot(error11,'LineWidth',1.5);hold on
plot(error1,'LineWidth',1.5)
legend('Haar','Linear')
xlabel('iteration')
ylabel('error')
title('WF PDHG')

subplot(235);
plot(error21,'LineWidth',1.5);hold on
plot(error2,'LineWidth',1.5)
legend('Haar','Linear')
xlabel('iteration')
ylabel('error')
title('WF ADMM')

subplot(236);
plot(error51,'LineWidth',1.5);hold on
plot(error5,'LineWidth',1.5)
legend('Haar','Linear')
xlabel('iteration')
ylabel('error')
title('WF Adaptive PDHG')



%%

[~,FDr1] = max(u1,[],2);
FDr1 = FDr1-1;
[uc,FDr2] = max(u2,[],2);
FDr2 = FDr2-1;

subplot(131);
scatter(D(1,:),D(2,:),5,FD);title('Ground Truth');axis square;
subplot(132);
scatter(D(1,:),D(2,:),5,FDr1);title('PDHG');axis square;
subplot(133);
scatter(D(1,:),D(2,:),5,FDr2);title('ADMM');axis square;
%}