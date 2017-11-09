function [unew, energy,residual,error]=TV_PDHGm_ClassK(FD0,Iset,u00,lambda,tol,G,maxit,adap_para,FD_ref)
% G is weight  matrix
% min lambda \sum_i|u_i|_WTV
% s.t u_i=FD0(:,i), sum (u_i)=1 i_i\geq 0

%%%%%%%%%%%%%5
% u=zeros(size(forig));
tic;

[M,K]=size(u00);
% find the index of nonzero G
%[r,c]=find(G);

Du=cell(K,1);
p0=cell(K,1);

%Isetc=setdiff(1:M, Iset);

% inititalize d and b, and compute normg
normtvg=0;
norm0=norm(u00,1);
% compute the initial  gradient and the norm
for k=1:K
    Du{k}=GraphGradientOperator(G,u00(:,k));
    normtvg=normtvg+sum(sum(G.*Du{k}));
    p0{k}=zeros(size(Du{k}));
end

%p0=Du;

energy=zeros(maxit,1);
residual=zeros(maxit,1);
error=zeros(maxit,1);

uold=u00;
unew=u00;
it=1;
stop=0;


theta=1;
%% parameters can be tuned
gamma=0.01;
sigma=sqrt(1/10)*10;
tau=sqrt(1/10)/10;

while (it<=maxit&&stop==0)
    %% update p
    ubar=unew+theta*(unew-uold);
    for k=1:K
        Du{k} = GraphGradientOperator(G,ubar(:,k));
        p0{k} = p0{k}+sigma*Du{k}; % old p is not stored due to the large memory
        % projection onto C_W ball
        % p0{k} = proj_W(p0{k},G*lambda);
        p0{k} = proj_W(p0{k},ones(size(G))*lambda);
    end
    %% update u
    uold = unew;
    for k=1:K
        unew(:,k) = uold(:,k)-tau*GraphGradientOperatorTranspose(G,p0{k});
        unew(Iset,k) = FD0(Iset,k);
    end
    % projection onto l1 ball
    unew = projl1p_1D(unew,1);
    
    % update  parameter: optional
    if (adap_para==1)
        theta = 1/sqrt(1+2*gamma*tau);
        tau = theta*tau;
        sigma = sigma/theta;
    end
    
    % compute the energy and residual
    residual(it) = norm(unew-uold,1)/norm0;
    for k=1:K
        energy(it) = energy(it)+sum(sum(G.*abs(Du{k})));
    end
    % compute the error if FD_ref is given
    [~,FDr]=max(unew,[],2);
    FDr=FDr-1;
    FDr(Iset)=FD_ref(Iset);
    %error(it)=sum(abs(FDr-FD_ref))/(M-length(Iset))*100;  %for 2 classes
    c=FDr==FD_ref;
    error(it)=100*(M-sum(c))/(M-length(Iset));
    if (residual<tol) 
        stop=1; 
    end
    if mod(it,20)==0
        Tm=toc;
        display( num2str(sum(abs(unew(:)-uold(:)))/length(unew)) );
        display(['Step = ' num2str(it) '; Residual = ' num2str(residual(it)) '; Energy = ' num2str(energy(it)) '; Accuracy = ' num2str(100-error(it)) '%; Time Elapsed = ' num2str(Tm)]);
    end
    it=it+1;
end

% display('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
% display('Program Finished.')
% display(['Step = ' num2str(nstep) '; Residual = ' num2str(residual) '; Error = ' num2str(error) '%; Time Elapsed = ' num2str(Tm)]);