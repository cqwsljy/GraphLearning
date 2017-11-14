function [u, energy,residual,error]=SplitBregGraphClassK2(FD0,Iset,u00,mu,lambda,dd,tol,W,WT,maxit,FD_ref)
% using K labelling functions
% min lambda \sum_i|Wu_i| .s.t u_i=FD0 on Ise,sum (u_i)=1 i_i\geq 0
tic;

[M,K] = size(u00);  %
d = cell(K,1);
b = cell(K,1);
Wu = cell(K,1);
deltab = cell(K,1);
mu2 = 1;
Isetc=setdiff(1:M, Iset);

% inititalize d and b, and compute normg
normg = 0;
for k=1:K
    d{k} = W(u00(:,k));
    b{k} = d{k};
    normg = normg+CoeffOperGraph('norm2',d{k});
    [dtmp,~] = find(FD0(:,k));
    eval(['Iset',num2str(k),'=dtmp;']);
    eval(['du',num2str(k),'=1./(dd(dtmp)+mu2);']);
end
[r Level] = size(d{1});
Thresh=cell(r,Level);
w=cell(r,Level);
for l=1:Level
    for j=1:r
        if (j == Level && j == 1)
            w{j,l} = lambda*4^(-l+1)*dd;
            Thresh{j,l} = w{j,l}/mu;
        else
            w{j,l} = lambda*4^(-l+1)*dd;
            Thresh{j,l}=w{j,l}/mu;
        end
    end
end
energy = zeros(maxit,1);
residual = zeros(maxit,1);
error = zeros(maxit,1);
u = u00;


disp(['Initial is ',num2str(100*length(Iset)/M),'%'])
for nstep=1:maxit

    % update u
    for k=1:K
        WTdb =  WT(CoeffOperGraph('-',d{k},b{k}));
        u(Isetc,k) = WTdb(Isetc);
        eval(['Isetmp=','Iset',num2str(k),';']);
        eval(['u(Isetmp,k) = du',num2str(k),'.*( FD0(Isetmp,k) + mu2*WTdb(Isetmp) );']);
    end
    % projection onto l1 ball
    u = projl1p_1D(u,1);

    % update d and b
    for k=1:K
        Wu{k} = W(u(:,k));%
        d{k} = CoeffOperGraph('s',CoeffOperGraph('+',Wu{k},b{k}),Thresh);
        deltab{k} = CoeffOperGraph('-',Wu{k},d{k});
        b{k} = CoeffOperGraph('+',b{k},deltab{k});

        % compute the energy and residual
        residual(nstep) = residual(nstep)+CoeffOperGraph('norm2',deltab{k})/normg;
        energy(nstep) = energy(nstep)+CoeffOperGraph('wnorm1',Wu{k},w);
    end

    % compute the error if FD_ref is given
    [~,FDr] = max(u,[],2);%
    FDr = FDr-1;FDr(Iset) = FD_ref(Iset);
    c = FDr == FD_ref;
    error(nstep) = 100*(M-sum(c))/(M - length(Iset));

    if residual(nstep)<tol
        break;
    end
    Tm=toc;
    if mod(nstep,20)==0
        display(['Step = ' num2str(nstep) '; Residual = ' num2str(residual(nstep)) '; Energy = ' num2str(energy(nstep)) '; Accuracy = ' num2str(100-error(nstep)) '%; Time Elapsed = ' num2str(Tm)]);
    end
end
toc