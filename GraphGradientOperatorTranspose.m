function DTp=GraphGradientOperatorTranspose(G,p)
% A: adacent matrix (diagonal is 0)  % normally p is the same as A 
% p: same size as A p(i,j) denote  u(i)-u(j) 

A = spones(G);
N = size(A,1);
Ap = A.*p; 
DTp = sum(Ap,1)'-sum(Ap,2);