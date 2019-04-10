function [x,lambda] = recycleRangeSpace(n,ubar,d0)

[H,g,A,b] = recycleConstruct(n,ubar,d0);
KKT = recycleKKT(n,ubar,d0);

L = chol(H,'lower');

v = (L*L')\g;

H_a = A'*(H)*A;
L_a = chol(H_a,'lower');
lambda = (L_a*L_a')\(b+A'*v);
x = H\(A*lambda-g);


end