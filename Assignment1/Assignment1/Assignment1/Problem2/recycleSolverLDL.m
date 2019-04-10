function [x,lambda] = recycleSolverLDL(n,ubar,d0)


[KKT,rhs] = recycleKKT(n,ubar,d0);

[L,D,p] = ldl(KKT,'lower','vector');  %Use LDL factorizatin because Matrix is symmetric and indefinite
sol(p) = L'\( D \(L\rhs(p)));


x = sol(1:n+1);
lambda = sol(n+2:end);

end