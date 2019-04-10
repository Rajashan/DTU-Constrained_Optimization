function [x,lambda,P] = recycleSparseLDLSolver(n,ubar,d0)

[KKT,rhs] = (recycleKKT(n,ubar,d0));

sp_KKT = sparse(KKT);

[L,D,P,Q] = ldl(sp_KKT,'lower','vector');

sol(P) = L'\(D\(L\rhs(P)));

x = sol(1:n+1);
lambda = sol(n+2:end);

end

