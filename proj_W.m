function p=proj_W(p0,W)
% projection onto l_infty<W_{ij} ball for each K

% % projection step
% tic
% p=p0;
% ind=find(abs(p0)>W+eps);
% % ind=find(abs(p0)>W);
% p(ind)=sign(p0(ind)).*W(ind);    
% toc

N = size(W,1);
p = p0(:);
W = W(:);

[r,~] = find(W);
ind = find( abs(p(r)) > W(r) );
p(r(ind)) = sign( p(r(ind)) ).*W(r(ind));
% t = full(sign( p(r(ind)) ).*W(r(ind)));
% for i = 1:length(ind)
%     j = ind(i);
%     p(r(j)) = t(i); 
%     t(i);
% end
p = reshape(p,N,N);