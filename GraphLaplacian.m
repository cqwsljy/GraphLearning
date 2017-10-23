function [L]=GraphLaplacian(D,h,K,graph_method)
if (nargin<2)
    h=10;
end

if (nargin<3)
    K=5;% K nearest neighboring vertices;
end

if nargin<4
    graph_method='knn';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compue adjacency matrix A
[~, M]=size(D);
A=zeros(M,M);
switch graph_method
    case 'knn'
        tree=vl_kdtreebuild(D);
        for k=1:M
            [id, dists] = vl_kdtreequery(tree, D, D(:,k),'NUMNEIGHBORS',K);
            A(k,id)=exp(-dists/h);
            %C(k,id)=1;
        end
    case 'full'  %full connected
        %%{
        clear C;
        for k=1:M
            Av=exp(-sum((D(:,k)*ones(1,M)-D).^2,1)/h);
            A(k,:)=Av;
        end
end
d=A*ones(M,1);
% L=eye(M,M)-diag(1./b)*A;
L=diag(d)-A;
end