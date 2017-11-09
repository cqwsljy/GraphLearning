function Du=GraphGradientOperator(A,u)
% u: n*1 vector
% A: adacent matrix (diagonal is 0) 
% compute the gradient A_ij (u(j)-u(i))
[r,c] = find(A);
N = size(A,1);

s = u(c)-u(r);
Du = sparse(r,c,s,N,N);
% Du=u(r)
% for i=1:N
%     UUj(c(r==i),:)=u(c(r==i))';
%     UUi(r(c==i),i)=u(r(c==i));
% end
% 
% % UUj=repmat(u',[N,1]); % replicate u matrix N rows
% % UUi=repmat(u,[1 N]); % replicate u matrix N columns
% 
% Du=sparse(A.*abs(UUj-UUi));