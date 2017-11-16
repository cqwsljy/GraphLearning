function DTp = GraphGradientOperatorTranspose(G,p)
% A: adacent matrix (diagonal is 0)  % normally p is the same as A 
% p: same size as A p(i,j) denote  u(i)-u(j) 
% Transpose of gradient is divergence,
% div(\phi)(x) = 1/(2d^r) \sum_(y \in v) w^q(\phi(x,y) - \phi(y,x))
% A = spones(G);
% N = size(A,1);
% Ap = A.*p; 
% DTp = sum(Ap,1)'- sum(Ap,2);

% Modified by Jiayong Liu,2017.11.14
Gp = G.*p; 
DTp = sum(Gp,1)'- sum(Gp,2);