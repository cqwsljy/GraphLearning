function [unew, energy,residual,errors]=TV_PDHGm_K(FD0,Iset,u00,lambda,tol,G,maxit,adap_para,FD_ref)
% G is weight  matrix
% min lambda \sum_i|u_i|_WTV
% s.t u_i=FD0(:,i), sum (u_i)=1 i_i\geq 0

%%%%%%%%%%%%%5
% u=zeros(size(forig));
tic;

[M,K] = size(u00);
% find the index of nonzero G
%[r,c]=find(G);

Du = cell(K,1);
p0 = cell(K,1);

%Isetc=setdiff(1:M, Iset);

% inititalize d and b, and compute normg
%normtvg=0;
norm0=norm(u00,1);
% compute the initial  gradient and the norm
for k=1:K
    Du{k} = GraphGradientOperator(G,u00(:,k));
    % normtvg=normtvg+sum(sum(G.*Du{k}));
    p0{k} =  Du{k};
end

%p0=Du;

energy = zeros(maxit,1);
residual = zeros(maxit,1);
errors = zeros(maxit,1);

% uold = zeros(M,K);
uold = rand(size(u00));
unew = uold;
%uold = u00;
%unew = u00;

it = 1;
stop = 0;


theta = 1;
%% parameters can be tuned
gamma = 0.01;

% sigma = sqrt(1/10)/100;
% tau = sqrt(1/10)*100;

% sigma = 0.02;% setpsize for dual variable
% tau = 10;  % setpsize for primal variable

sigma = 0.001;% setpsize for dual variable
tau = 300;  % setpsize for primal variable


while (it<=maxit && stop== 0)
    %% update p
    ubar = unew + theta*(unew - uold);
    for k=1:K
        Du{k} = GraphGradientOperator(G,ubar(:,k)); 
        p0{k} = p0{k} + sigma*Du{k}; % old p is not stored due to the large memory
        % projection onto C_W ball
        p0{k} = proj_W(p0{k},G*lambda);
        % p0{k} = proj_W(p0{k},ones(size(G))*lambda);
        % p0{k} = p0{K}.*(abs(p0{k}) <= 1) + sign(p0{K}) .* (abs(p0{k}) > 1);
    end
    
    %% update u
    uold = unew;
    for k=1:K
        unew(:,k) = uold(:,k) - tau*GraphGradientOperatorTranspose(G,p0{k});
        %unew(Iset(:,k),k) = FD0(Iset(:,k),k);
    end

    for k = 1:K
        unew(Iset(:,k),:) = 0; 
    end
    for k = 1:K
        unew(Iset(:,k),k) = FD0(Iset(:,k),k);
    end

    % projection onto l1 ball
    unew = projl1p_1D(unew,1);
    
    % update  parameter: optional
    %{
    if (adap_para==1 && it > 100)
        theta = 1/sqrt(1+2*gamma*tau);
        tau = theta*tau;
        sigma = sigma/theta;
    end
    %}
    
    % compute the energy and residual
    residual(it) = norm(unew-uold,1)/norm0;
    for k=1:K
        Du{k} = GraphGradientOperator(G,unew(:,k));
        energy(it) = energy(it) + sum(sum(G.*abs(Du{k})));
    end
    
    % compute the error if FD_ref is given
    [~,FDr]=max(unew,[],2);
    FDr = FDr-1;
    FDr(Iset(:)) = FD_ref(Iset(:));
    c = FDr == FD_ref;
    errors(it)=100*(M-sum(c))/(M-length(Iset(:)));
    if (residual(it) < tol) 
        errors = errors(1:it);
        residual = residual(1:it);
        energy = energy(1:it);
        stop=1; 
    end
    if mod(it,2)==0
        Tm=toc;
        display( num2str(sum(abs(unew(:)-uold(:)))/length(unew)) );
        display(['Step = ' num2str(it) '; Residual = ' num2str(residual(it)) '; Energy = ' num2str(energy(it)) '; Accuracy = ' num2str(100-errors(it)) '%; Time Elapsed = ' num2str(Tm)]);
    end
    it=it+1;
end