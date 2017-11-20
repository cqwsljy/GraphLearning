function [y] = soft(x,alpha)
y = sign(x).*max(abs(x)-alpha,0);
end