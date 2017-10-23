function [u, energy,residual,error]=SplitBregGraphClassKexact(FD0,Iset,u00,mu,mu1,lambda,dd,tol,W,WT,maxit,FD_ref)
% using K labelling functions
% min lambda \sum_i|Wu_i| .s.t u_i=FD0 on Ise,sum (u_i)=1 i_i\geq 0
tic;
% u=zeros(size(FD_ref));
[M,K] = size(u00);  %
d = cell(K,1);
b = cell(K,1);
b1 = zeros(length(Iset),K);
Wu = cell(K,1);
deltab = cell(K,1);

Isetc=setdiff(1:M, Iset);
% inititalize d and b, and compute normg
normg = 0;
for k=1:K
    Wuk = W(u00(:,k)); 
    d{k} = Wuk; 
    b{k} = Wuk;
    normg = normg+CoeffOperGraph('norm2',Wuk);
end
[r, Level] = size(d{1});
% Thresh=cell{r,Level};
% w=cell{r,Level};
for l=1:Level
    for j=1:r
        w{j,l} = lambda*4^(-l+1)*dd;
        Thresh{j,l} = w{j,l}/mu;
    end
end
energy = zeros(maxit,1);
residual = zeros(maxit,1);
error = zeros(maxit,1);
u = u00;
ut = u;
disp(['Initial is ',num2str(100*length(Iset)/M),'%'])
for nstep=1:maxit
    %
    for k=1:K
        WTdb =  WT(CoeffOperGraph('-',d{k},b{k}));
        u(Isetc,k) = WTdb(Isetc);
        %u(Iset,k) = FD0(Iset,k);
        u(Iset,k) = 1./(mu1 + mu) .* (mu1*(FD0(Iset,k)-b1(:,k))+mu*WTdb(Iset));
    end
    %u(u<0)=0;u(u>1)=1;
    % projection onto l1 ball
    u = projl1p_1D(u,1);
    % update d and b
    for k=1:K
        Wu{k} = W(u(:,k));%
        % update d
        d{k} = CoeffOperGraph('s',CoeffOperGraph('+',Wu{k},b{k}),Thresh);
        % update b
        deltab{k} = CoeffOperGraph('-',Wu{k},d{k});
        b{k} = CoeffOperGraph('+',b{k},deltab{k});
        % update b1
        b1(:,k) = b1(:,k) + u(Iset,k) - FD0(Iset,k);
        % compute the energy and residual
        % residual(nstep) = residual(nstep)+CoeffOperGraph('norm2',deltab{k})/normg;
        energy(nstep) = energy(nstep)+CoeffOperGraph('wnorm1',Wu{k},w);
    end
    
    %residual(nstep) = sqrt( sum( sum((ut-u).^2,2) ) )/normg;
    if nstep > 1
        residual(nstep) = abs(energy(nstep) - energy(nstep -1));
    else
        residual(nstep) = 100;
    end
    ut = u;
    % energy(nstep) = energy(nstep)+sum(sum(abs(u(Iset,:)-FD0(Iset,:))));
    % compute the error if FD_ref is given
    [~,FDr] = max(u,[],2);%
    FDr = FDr-1;FDr(Iset) = FD_ref(Iset);
    c = FDr == FD_ref;
    %error(nstep) = 100*(M-sum(c))/(M-length(Iset));
    error(nstep) = 100*(M-sum(c))/(M - length(Iset));
    %error(nstep) = 100-sum(length(c(c==1))-length(Iset))/(M-length(Iset))*100;
    
    %disp(['Initial error is ',num2str(error),'%'])
    if residual(nstep)<tol
        break;
    end
    Tm=toc;
    if mod(nstep,200)==0
        display(['Step = ' num2str(nstep) '; Residual = ' num2str(residual(nstep)) '; Energy = ' num2str(energy(nstep)) '; Accuracy = ' num2str(100-error(nstep)) '%; Time Elapsed = ' num2str(Tm)]);
    end
end
toc
% display('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
% display('Program Finished.')
% display(['Step = ' num2str(nstep) '; Residual = ' num2str(residual) '; Error = ' num2str(error) '%; Time Elapsed = ' num2str(Tm)]);
