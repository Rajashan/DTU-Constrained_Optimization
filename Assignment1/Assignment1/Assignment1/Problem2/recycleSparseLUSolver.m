function [x,lambda,P] = recycleSparseLUSolver(n,ubar,d0)

[KKT,rhs] = (recycleKKT(n,ubar,d0));

sp_KKT = sparse(KKT);

[L,U,P] = lu(sp_KKT);

sol = U\(L\(P*rhs));

x = sol(1:n+1);
lambda = sol(n+2:end);

end

