function [u, nstep]=SplitBregGraphMultiClass(FD0,Iset,K,mu,lambda,tol,W,WT,maxit,forig)
% use one label function for multiclass
tic;
u=zeros(size(forig));
M=length(u);
d=W(u);b=W(u);
normg=CoeffOperGraph('norm2',W(forig));
[r Level]=size(d);
for l=1:Level
    for j=1:r
        Thresh{j,l}=lambda/mu*4^(-l+1);
    end
end
Isetc=setdiff(1:M, Iset);
for nstep=1:maxit
    WTdb=WT(CoeffOperGraph('-',d,b));
    u(Isetc)=WTdb(Isetc);
    u(Iset)=1/(1+mu)*(FD0+mu*WTdb(Iset));
%     u(u<0)=0;u(u>3)=3;
    Wu=W(u);
    d=CoeffOperGraph('s',CoeffOperGraph('+',Wu,b),Thresh);
%     d=CoeffOperGraph('s_band',CoeffOperGraph('+',Wu,b),Thresh,[3:r]);
    deltab=CoeffOperGraph('-',Wu,d);
    b=CoeffOperGraph('+',b,deltab);
    residual=CoeffOperGraph('norm2',deltab)/normg;
    ut=real(u);ut(Iset)=FD0;
	% projection onto the region
	ut(ut>K-1.5)=K-1;
	for i=1:K-1
     ut((ut<i+0.5)&(ut>=i-0.5))=i;
	 end
	ut(ut<0.5)=0;
    error=sum(double(ut~=forig))/(M-length(Iset))*100;
    if residual<tol
        break;
    end
    Tm=toc;
    if mod(nstep,200)==0
        display(['Step = ' num2str(nstep) '; Residual = ' num2str(residual) '; Error = ' num2str(error) '%; Time Elapsed = ' num2str(Tm)]);
        hist(u,100);xlim([0 3]);drawnow;
    end
end
display('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
display('Program Finished.')
display(['Step = ' num2str(nstep) '; Residual = ' num2str(residual) '; Error = ' num2str(error) '%; Time Elapsed = ' num2str(Tm)]);