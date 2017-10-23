function [u, energy,residual]=SplitBregGraphCluster_Potts(FD0,u00,mu,lambda,dd,tol,W,WT,maxit)
% min lambda \sum_i|Wu_i|+ sum_i <u_i, FD0>
tic;
% u=zeros(size(forig));
[M,K]=size(FD0);
d=cell(K,1);
b=cell(K,1);
Wu=cell(K,1);
deltab=cell(K,1);

% inititalize d and b, and compute normg
normg=0;
for k=1:K
d{k}=W(u00(:,k));
b{k}=d{k};
normg=normg+CoeffOperGraph('norm2',d{k});
end

[r Level]=size(d{1});
for l=1:Level
    for j=1:r
             w{j,l}=lambda*4^(-l+1)*dd;
            Thresh{j,l}=w{j,l}/mu;

end
end

WTdb=zeros(M,K);
energy=zeros(maxit,1);
residual=zeros(maxit,1);
u=u00;
for nstep=1:maxit
     % 
     for k=1:K
    WTdb=WT(CoeffOperGraph('-',d{k},b{k})); 
    u(:,k)=(-FD0(:,k)+mu*WTdb)/mu;  
     end
       %u(u<0)=0;u(u>1)=1;
       % projection onto l1 ball
       u = projl1p_1D(u,1);
     % update d and b
     for k=1:K
    Wu{k}=W(u(:,k));
     d{k}=CoeffOperGraph('s',CoeffOperGraph('+',Wu{k},b{k}),Thresh);
%     d=CoeffOperGraph('s_band',CoeffOperGraph('+',Wu,b),Thresh,[3:r]);
    deltab{k}=CoeffOperGraph('-',Wu{k},d{k});
    b{k}=CoeffOperGraph('+',b{k},deltab{k});
    % compute the energy and residual
    residual(nstep)=residual(nstep)+CoeffOperGraph('norm2',deltab{k})/normg;
    energy(nstep)=energy(nstep)+CoeffOperGraph('wnorm1',Wu{k},w);
     end
     energy(nstep)=energy(nstep)+sum(sum(u.*FD0));
     if residual<tol
        break;
    end
    Tm=toc;
    if mod(nstep,100)==0
        display(['Step = ' num2str(nstep) '; Residual = ' num2str(residual(nstep)) '; Energy = ' num2str(energy(nstep)) '; Time Elapsed = ' num2str(Tm)]);
    end
end
% display('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
% display('Program Finished.')
% display(['Step = ' num2str(nstep) '; Residual = ' num2str(residual) '; Error = ' num2str(error) '%; Time Elapsed = ' num2str(Tm)]);