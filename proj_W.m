function p=proj_W(p0,W)
% projection onto l_infty<W_{ij} ball for each K

% projection step
p=p0;
ind=find(abs(p0)>W+eps);
p(ind)=sign(p0(ind)).*W(ind);    