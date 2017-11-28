function unew = CGforTV(L,G,uold,d,b,FD0,Iset)
format long
Classk = size(uold,2);
n = size(L,1);
tol = 1e-10;
r = cell(Classk,1);
p = cell(Classk,1);
r2 = cell(Classk,1);
Du = cell(1,Classk);
for k = 1:Classk
    r{k} = GraphGradientOperatorTranspose(G,d{k} - b{k}) - L*uold(:,k);
    %r{k} = b{k} - L*uold(:,k);
    p{k} = r{k};  
end
unew = uold;
energy = zeros(Classk,1);
for i = 0:(n-1)
    uold = unew;
    for k = 1:Classk
        alpha = (r{k}'*r{k})/(p{k}'*L*p{k});
        unew(:,k) = uold(:,k) + alpha*p{k};
        r2{k} = r{k} - alpha*L*p{k};
        if ((norm(r2{k}) <= tol)||(i == n-1))
           continue;
        end
        beta = norm(r2{k})^2/norm(r{k})^2;
        p{k} = r2{k} + beta*p{k};
        r{k} = r2{k};
%         disp(num2str(norm(unew(:) - uold(:))/length(unew)));
%         disp(num2str(norm(r2{k})));
%        Du{k} = GraphGradientOperator(G,unew(:,k)); % gradient at u^{i+1}
%        energy(k) = sum(sum(G.*abs(Du{k})))
    end
    %%{
    for k = 1:Classk
        unew(Iset(:,k),:) = 0; 
    end
    for k = 1:Classk
        unew(Iset(:,k),k) = FD0(Iset(:,k),k);
    end
    %}

    flag = 1;

    for k = 1:Classk
        flag = flag * (norm(r2{k}) <= tol);
        %energy = unew(:,k)'*L'*L*unew(:,k) - 2*unew(:,k)'*L*GraphGradientOperatorTranspose(G,d{k} - b{k});
    end
    % disp(num2str(sum(abs(unew(:) - uold(:)))/length(uold)));
    %disp(num2str(sum(energy)));
%    if (mod(i,100) ==0)
%        disp(num2str(sum(abs(unew(:) - uold(:)))/length(uold)));
%    end

    unew = projl1p_1D(unew,1);

    if (flag  == 1)
        break
    end
end
disp(num2str(sum(abs(unew(:) - uold(:)))/length(uold)));