function [L,d,lambda_max]=GenerateGraph_fun(D,h,K,graph_method)
% generate the graph laplacian matrix L
% degree d;
% One Eigenvector corresponds to the largest eigenvalue
% Largest Eigenvalue
%revised 5/6/2016
%%%%%%%%%%%%%%%%% Synthetic Data %%%%%%%%%%%%%%%%%%%
% load('.\Synthetic\DataSyn_3moons.mat'); % h=10; K=5;
% load('.\Synthetic\DataSyn_2Circles.mat'); % h=10; K=5;
% load('.\Synthetic\DataSyn_2Circles_p1000.mat'); % h=10; K=5;
%load('.\Synthetic\DataSyn_4Circles.mat'); % h=10; K=5;
% load('.\Synthetic\DataSyn_Annulus.mat'); % h=10; K=5;
% load('.\Synthetic\DataSyn_Annulus_p1000.mat'); % h=10; K=5;

%%%%%%%%%%%%%%%%%%%% Real Data %%%%%%%%%%%%%%%%%%%%%
if (nargin<2)
    h=10;
end

if (nargin<3)
    K=5;% K nearest neighboring vertices;
end

if nargin<4
    graph_method='knn';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

[~, M]=size(D);
A=zeros(M,M);   %weighted matrix
C=zeros(M,M);   %Logic matrix to construct symetric Laplacian matrix
switch graph_method
    case 'knn'
        %original code
        %{
        for k=1:M
            Av=exp(-sum((D(:,k)*ones(1,M)-D).^2,1)/h);
            [~, I]=sort(Av,'descend');
            A(k,I(1:K+1))=Av(I(1:K+1));
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        %use kdtree1.2,knn
        D=D';
        tree = kdtree_build(D);
        for k=1:M
            [idxs1, ~] = kdtree_k_nearest_neighbors(tree,D(k,:),K);
            id=flipud(idxs1);
            %Av=flipud(dists);
            Av=exp(-sum((repmat(D(k,:),K,1)-D(id,:)).^2,2)/h);
            A(k,id)=Av;
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%{
        %use vlfead package,knn,symmetric
        tree=vl_kdtreebuild(D);
        for k=1:M
            [id, dists] = vl_kdtreequery(tree, D, D(:,k),'NUMNEIGHBORS',K);
            A(k,id)=exp(-dists.^2/h);
            %C(k,id)=1;
        end
        %C=C.*C';
        %A=A.*C;%symmerticlized
    case 'full'  %full connected
        %%{
        clear C;
        for k=1:M
            Av=exp(-sum((D(:,k)*ones(1,M)-D).^2,1)/h);
            A(k,:)=Av;
        end
        %}
    case 'ZM' 
        % Zelnik-Manor and Perona weight function
        % w(x,y) = exp(-d(x,y)^2/sqrt(tau(x)*tau(y)))
        % where sqrt(tau(x)) is the distance between x and its Kth closest neigthbor
        tree = vl_kdtreebuild(D);
        h2 = zeros(K,1);
        distnearK = zeros(size(D,2),1);
        for i = 1:M
        	[~, dists] = vl_kdtreequery(tree,D,D(:,i),'NUMNEIGHBORS',K);
        	distnearK(i) = dists(end);
        end
        for i = 1:M
            [id, dists] = vl_kdtreequery(tree,D,D(:,i),'NUMNEIGHBORS',K);
            dx = sqrt(sum((D(:,i) - D(:,id(end))).^2));
            for k = 1:K
                dy = distnearK(id(k));
                h2(k) = dx * dy;
            end
            A(i,id) = exp(-dists.^2./h2);
        end
    case 'ZM2' 
        % Zelnik-Manor and Perona weight function
        % w(x,y) = exp(-d(x,y)^2/sqrt(tau(x)*tau(y)))
        % where sqrt(tau(x)) is the distance between x and its Kth closest neigthbor
        tree = vl_kdtreebuild(D);
        h2 = zeros(K,1);
        distnearK = cell(size(D,2),1);
        idnearK = cell(size(D,2),1);

        for i = 1:M
            [id, dists] = vl_kdtreequery(tree,D,D(:,i),'NUMNEIGHBORS',K);
            distnearK{i} = dists;
            idnearK{i} = id;
        end
        for i = 1:M
            dists = distnearK{i};
            id = idnearK{i};
            dx = sqrt(sum((D(:,i) - D(:,id(end))).^2));
            for k = 1:K
                dists = distnearK{id(k)};
                dy = dists(end);
                h2(k) = dx * dy;
            end
            A(i,id) = exp(-dists.^2./h2);
        end
    case 'ZM3' 
        % Zelnik-Manor and Perona weight function
        % w(x,y) = exp(-d(x,y)^2/sqrt(tau(x)*tau(y)))
        % where sqrt(tau(x)) is the distance between x and its Kth closest neigthbor
        tree = vl_kdtreebuild(D);
        h2 = zeros(K,1);
        for i = 1:M
            [id, dists] = vl_kdtreequery(tree,D,D(:,i),'NUMNEIGHBORS',K);
            dx = sqrt(sum((D(:,i) - D(:,id(end))).^2));
            for k = 1:K
                [~, distsy] = vl_kdtreequery(tree,D,D(:,id(k)),'NUMNEIGHBORS',K);
                dy = distsy(end);
                h2(k) = dx * dy;
            end
            A(i,id) = exp(-dists.^2./h2);
        end
end
% symmetric
A = max(A,A');
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construct Laplacian matrix
d=A*ones(M,1);
% L=eye(M,M)-diag(1./b)*A;
L=diag(d)-A;
L=sparse(L);
lambda_max=eigs(L,1);%
%[Usm Lambda_sm]=eigs(L,2,'sm'); % error for some data, need to fix
%U1=Usm(:,2);
toc