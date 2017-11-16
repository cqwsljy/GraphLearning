function unew = CGforTV(L,G,uold,d,b,FD0,Iset)
Classk = size(uold,2);
n = size(L,1);
tol = 1e-4;
r = cell(Classk,1);
y = cell(Classk,1);
for k = 1:Classk
    r{k} = GraphGradientOperatorTranspose(G,d{k} - b{k}) - L*uold(:,k);
    y{k} = r{k};  
end

for i = 0:(n-1)
    unew = uold;
    for k = 1:Classk
        alpha = (r{k}'*r{k})/(y{k}'*L*y{k});
        unew(:,k) = uold(:,k)+alpha*y{k};
        r2 = GraphGradientOperatorTranspose(G,d{k} - b{k}) - L*uold(:,k); 
        if ((norm(r2) <= tol)||(k == n-1))
           break;
        end
        beta = norm(r2)^2/norm(r{k})^2;
        y{k} = r2 + beta*y{k};
        r{k} = r2;
        uold(Iset,k) = FD0(Iset,k);
%         disp(num2str(norm(unew(:) - uold(:))/length(unew)));
        disp(num2str(norm(r2)));
    end
end

% function x = cg(A,b)
% tol=1e-10;
% r = b + A*b;
% w = -r;
% z = A*w;
% s = w'*z;
% t = (r'*w)/s;
% x = -b + t*w;
% for k = 1:numel(b);
%     r = r - t*z;
%     if( norm(r) < tol )
%        return;
%     end
%     B = (r'*z)/s;
%     w = -r + B*w;
%     z = A*w;
%     s = w'*z;
%     t = (r'*w)/s;
%     x = x + t*w;
% end

% function x=Gongetidu2(A,b,x0,epsa)
% n=size(A,1);
% x=x0;
% r=b-A*x;
% d=r;
% for k=0:(n-1)
%     alpha=(r'*r)/(d'*A*d);
%     x=x+alpha*d;
%     r2=b-A*x; 
%     if ((norm(r2)<=epsa)||(k==n-1))
%        break;
%     end
%     beta=norm(r2)^2/norm(r)^2;
%     d=r2+beta*d;
%     r=r2;
% end
% % 