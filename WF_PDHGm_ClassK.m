function [unew, energy,residual,error] = WF_PDHGm_ClassK(FD0,Iset,u00,lambda,dd,tol,W,WT,maxit,adap_para,FD_ref)
% using K labelling functions
% min lambda \sum_i|Wu_i| .s.t u_i=FD0 on Iset; u in Simplex.
tic;
% u=zeros(size(forig));
[M,K]=size(u00);
d=cell(K,1);
Wu=cell(K,1);
% inititalize d and b, and compute normg
normg=0;
norm0=norm(u00,1);

for k=1:K
    d{k}=W(u00(:,k));
    normg=normg+CoeffOperGraph('norm2',d{k});
end

[r Level]=size(d{1});
w = cell(r,Level);
Thresh = w;
for l=1:Level
    for j=1:r
        if (j == 1 && l == Level)
            w{j,l} = 0*lambda*4^(-l+1)*dd;
            Thresh{j,l} = w{j,l};
        else
            w{j,l} = lambda*4^(-l+1)*dd;
            Thresh{j,l} = w{j,l};
        end
    end
end

energy=zeros(maxit,1);
residual=zeros(maxit,1);
error=zeros(maxit,1);
uold=zeros(M,K);
unew=uold;

theta=1;
% sigma, tau can be tuned
gamma = 0.01; % parameter for update theta
sigma = 0.03;% setpsize for dual variable
tau = 50;  % setpsize for primal variable

for nstep=1:maxit
    ubar=unew+theta*(unew-uold);
    % update d
    for k=1:K
        Wu{k}=W(ubar(:,k));
        d{k}=CoeffOperGraph('*+',d{k},Wu{k},1,sigma); % compute d=d+sigma*Wu;
        % d^{k+1} = (I + sigma \partial J^*) ^{-1} * ()
        d{k} = CoeffOperGraph('p',d{k},Thresh); % projection onto l infinity ball with Thresh
    end
    
    % update u
    uold = unew;
    for k=1:K
        unew(:,k)=uold(:,k)-tau*WT(d{k});
        unew(Iset,k)=FD0(Iset,k);
    end
    
    % projection onto l1 ball
    unew = projl1p_1D(unew,1);
    
    % update  parameter: optional
    %%{
    if (adap_para==1 && nstep > 300)
        theta = 1/sqrt(1+2*gamma*tau);
        tau = theta*tau;
        sigma = sigma/theta;
    end
    %}
    
    % Compute the enery and residual
    residual(nstep)=norm(unew-uold,1)/norm0;
    for k=1:K
        energy(nstep)=energy(nstep)+CoeffOperGraph('wnorm1',Wu{k},Thresh);
    end
    
    % compute the error if FD_ref is given
    [~,FDr] = max(unew,[],2);
    FDr = FDr-1;FDr(Iset) = FD_ref(Iset);
    c = FDr == FD_ref;
    error(nstep) = 100*(M-sum(c))/(M - length(Iset));%sum(c) contains length(Iset),so not minus length(Iset) on numerator
    if residual<tol
        break;
    end
    if mod(nstep,20)==0
        Tm=toc;
        display(['Step = ' num2str(nstep) '; Residual = ' num2str(residual(nstep)) '; Energy = ' num2str(energy(nstep)) '; Accuracy = ' num2str(100-error(nstep)) '%; Time Elapsed = ' num2str(Tm)]);
    end
end
% display('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
% display('Program Finished.')
% display(['Step = ' num2str(nstep) '; Residual = ' num2str(residual) '; Error = ' num2str(error) '%; Time Elapsed = ' num2str(Tm)]);