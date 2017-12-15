function [unew, energy,residual,errors] = WF_PDHGm_ClassK(FD0,Iset,u00,lambda,dd,mu,tol,W,WT,maxit,adap_para,FD_ref)
% using K labelling functions
% min lambda \sum_i|Wu_i| .s.t u_i=FD0 on Iset; u in Simplex.
tic;
% u=zeros(size(forig));
[M,K]=size(u00);
d=cell(K,1);
Wu=cell(K,1);
Wue = cell(K,1);
% inititalize d and b, and compute normg
normg = 0;
norm0 = norm(u00,1);

for k=1:K
    d{k} = W(u00(:,k));
    normg = normg+CoeffOperGraph('norm2',d{k});
end


[r Level]=size(d{1});
w = cell(r,Level);
Thresh = w;
for l=1:Level
    for j=1:r
        if (j == 1 && l == Level)
            w{j,l} = lambda*4^(-l+1)*dd;
            Thresh{j,l} = w{j,l}/mu;
        else
            w{j,l} = lambda*4^(-l+1)*dd;
            Thresh{j,l} = w{j,l}/mu;
        end
    end
end


energy = zeros(maxit,1);
residual = zeros(maxit,1);
errors = zeros(maxit,1);
uold = randn(size(u00));
% uold = zeros(M,K);
% uold = u00;
unew = uold;

theta = 1;
% sigma, tau can be tuned
gamma = 0.01; % parameter for update theta
sigma = 0.008;% setpsize for dual variable
tau = 30;  % setpsize for primal variable

disp(['Initial is ',num2str(100*length(Iset(:))/M),'%'])
for nstep=1:maxit
    ubar = unew + theta*(unew-uold);
    uold = unew;
    % update d and u
    for k=1:K
        % update d
        Wu{k}=W(ubar(:,k));
        d{k} = CoeffOperGraph('*+',d{k},Wu{k},1,sigma); % compute d=d+sigma*Wu;
        % d^{k+1} = (I + sigma \partial J^*) ^{-1} * ()
        d{k} = CoeffOperGraph('p',d{k},Thresh); % projection onto l infinity ball with Thresh

        % update u
        unew(:,k)=uold(:,k)-tau*WT(d{k});
        % unew(Iset,k)=FD0(Iset,k);
        % unew(Iset(:,k),k)=FD0(Iset(:,k),k);
        
    end
    % projection onto l1 ball
    %%{
    % unew = projl1p_1D(unew,1);

    unew = projl1p_1D(unew,1);
    for k = 1:K
        unew(Iset(:,k),:) = 0; 
    end
    for k = 1:K
        unew(Iset(:,k),k) = FD0(Iset(:,k),k);
        Wue{k} = W(unew(:,k));
    end
    %}

    % update  parameter: optional
    %{
    if (adap_para==1 && nstep > 0)
        theta = 1/sqrt(1+2*gamma*tau);
        tau = theta*tau;
        sigma = sigma/theta;
    end
    %}

    
    % Compute the enery and residual
    residual(nstep)=norm(unew-uold,1)/norm0;
    for k=1:K
        Wue{k}=W(unew(:,k));
        energy(nstep)=energy(nstep) + CoeffOperGraph('wnorm1',Wue{k},w);
    end
    
    % compute the errors if FD_ref is given
    [~,FDr] = max(unew,[],2);
    FDr = FDr-1;FDr(Iset(:)) = FD_ref(Iset(:));
    c = FDr == FD_ref;
    errors(nstep) = 100*(M-sum(c))/(M - length(Iset(:)));%sum(c) contains length(Iset),so not minus length(Iset) on numerator
    if residual(nstep)<tol
        errors = errors(1:nstep)
        residual = residual(1:nstep)
        energy = energy(1:nstep)
        break;
    end
    if mod(nstep,1)==0
        Tm=toc;
        display(['Step = ' num2str(nstep) '; Residual = ' num2str(residual(nstep)) '; Energy = ' num2str(energy(nstep)) '; Accuracy = ' num2str(100-errors(nstep)) '%; Time Elapsed = ' num2str(Tm)]);
    end
end
% display('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
% display('Program Finished.')
% display(['Step = ' num2str(nstep) '; Residual = ' num2str(residual) '; errors = ' num2str(errors) '%; Time Elapsed = ' num2str(Tm)]);