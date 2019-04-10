function [x,lambda] = recycleSolverLU(n,ubar,d0)

[KKT,rhs] = recycleKKT(n,ubar,d0);


[L,U,p] = lu(KKT,'vector');

sol(p) = U\(L\rhs(p));

x = sol(1:n+1)';
lambda = sol(n+2:end)';


end





