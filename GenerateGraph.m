clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Classification %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% Synthetic Data %%%%%%%%%%%%%%%%%%%

% load('.\Synthetic\DataSyn_2Circles.mat'); % h=10; K=5;
% load('.\Synthetic\DataSyn_2Circles_p1000.mat'); % h=10; K=5;
load('.\Synthetic\DataSyn_4Circles.mat'); % h=10; K=5;
% load('.\Synthetic\DataSyn_Annulus.mat'); % h=10; K=5;
% load('.\Synthetic\DataSyn_Annulus_p1000.mat'); % h=10; K=5;

%%%%%%%%%%%%%%%%%%%% Real Data %%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

[tt M]=size(D);

h=10;
K=5;% K nearest neighboring vertices;
A=zeros(M,M);
for k=1:M
    Av=exp(-sum((D(:,k)*ones(1,M)-D).^2,1)/h);
    [tp I]=sort(Av,'descend');
    A(k,I(1:K+1))=Av(I(1:K+1));
end
d=A*ones(M,1);
% L=eye(M,M)-diag(1./b)*A;
L=diag(d)-A;
L=sparse(L);
lambda_max=eigs(L,1);
[Usm Lambda_sm]=eigs(L,2,'sm');%Smallest magnitude. Same as sigma = 0. 
							   %If A is a function, Afun must return Y = A\x.

U1=Usm(:,2);%第二小的特征值对应的特征向量

toc



