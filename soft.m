function [x] = soft(x,alpha,flag)
	if nargin > 2
		[M,N] = size(x);
		x = x(:);
		alpha = alpha(:);
		indexx = find(x~=0);
		indexy = find(alpha~=0);
		index = union(indexx,indexy);
		z1 = x(index);
		z2 = alpha(index);
		z = sign(z1).*max(abs(z1)-z2,0);
		x(index) = z;
        x = reshap(x,M,N);
	else
		x = sign(x).*max(abs(x)-alpha,0);
	end
end