function [unew, energy,residual,error]=TV_SplitBregClassK(FD0,Iset,u00,tol,G,maxit,mu,FD_ref)
% G is weight  matrix
% min lambda \sum_i|u_i|_WTV
% s.t u_i=FD0(:,i), sum (u_i)=1 i_i\geq 0

%%%%%%%%%%%%%5
% u=zeros(size(forig));
tic;

[M,K]=size(u00);
% find the index of nonzero G
%[r,c]=find(G);

Du = cell(K,1);
d = cell(K,1);
b = cell(K,1);
deltab = cell(K,1);
%Isetc=setdiff(1:M, Iset);

% inititalize d and b, and compute normg
%normtvg = 0;
norm0 = norm(u00,1);
uold = zeros(size(u00));
unew = zeros(size(u00));
% compute the initial  gradient and the norm
for k = 1:K
    Du{k} = GraphGradientOperator(G,uold(:,k));
    %normtvg = normtvg+sum(sum(G.*Du{k}));
    d{k} = rand(size(Du{k}));
    b{k} = rand(size(Du{k}));
end

%p0=Du;

energy=zeros(maxit,1);
residual=zeros(maxit,1);
error=zeros(maxit,1);


it=1;
stop=0;


%% parameters can be tuned
delta = 0.0071; % stepsize for least square in update u
dd = sum(G);
D = diag(dd);
L = D - G;
A = L'*L;
while (it<=maxit&&stop==0)
    %% update u
    flag = 1;
    while (flag)
        uold = unew;
        for k=1:K
            unew(:,k) = uold(:,k) - delta*( GraphGradientOperatorTranspose(G, GraphGradientOperator(G,uold(:,k)) - d{k} + b{k}) );
        end
%         for k = 1:K
%             unew(Iset(:,k),:) = 0; 
%         end
%         for k = 1:K
%             unew(Iset(:,k),k) = FD0(Iset(:,k),k);
%         end
%         unew = projl1p_1D(unew,1);
        residualU = norm(unew(:)-uold(:))/size(unew,1);
        disp(num2str(residualU));
        flag = residualU > 1e-7;
    end

%     uold = unew;
%     unew = CGforTV(L,G,uold,d,b,FD0,Iset);

    % projection onto l1 ball
    unew = projl1p_1D(unew,1);
%     uold = unew;
    %% update d
    
    for k=1:K
        Du{k} = GraphGradientOperator(G,unew(:,k)); % gradient at u^{i+1}
        d{k} = soft(Du{k} + b{k},G/mu);
        % tmp = Du{k} + b{k};
        %c1 = (tmp > G/mu);
        %c2 = (tmp <- G/mu);
        % c = (tmp <= G/mu) & (tmp >=- G/mu);
        %d{k}(c) = 0;
        %d{k}(c1) = G.*(tmp-1)/mu;
        %d{k}(c2) = G.*(tmp + 1)/mu; 
        %d{k} = (G.*(tmp-1)./mu).*(tmp > G/mu) + (G.*(tmp + 1)/mu).*(tmp < G/mu);
        %d{k}(c) = 0;
    end

    %% update b
    for k = 1:K
        deltab{k} = Du{k} - d{k};
        b{k} = b{k} + deltab{k};
    end

    
    % compute the energy and residual
    residual(it) = norm(unew-uold,1)/norm0;
    for k=1:K
        energy(it) = energy(it) + sum(sum(G.*abs(Du{k})));
        residual(it) = residual(it)+sum(sum(deltab{k}.^2))/size(unew,1);
    end

    % compute the error if FD_ref is given
    [~,FDr] = max(unew,[],2);
    FDr = FDr-1;
    FDr(Iset(:)) = FD_ref(Iset(:));
    c = FDr == FD_ref;
    error(it)=100*(M-sum(c))/(M-length(Iset(:)));
    if (residual<tol) 
        stop=1; 
    end
    if mod(it,2)==0
        Tm = toc;
        % display( num2str(sum(abs(unew(:)-uold(:)))/length(unew)) );
        display(['Step = ' num2str(it) '; Residual = ' num2str(residual(it)) '; Energy = ' num2str(energy(it)) '; Accuracy = ' num2str(100-error(it)) '%; Time Elapsed = ' num2str(Tm)]);
    end
    it = it+1;
end