function [KKT,rhs] = recycleKKT(n,ubar,d0)

[H,g,A,b] = recycleConstruct(n,ubar,d0);

KKT = [H, -A; -A',zeros(n,n)];
rhs = -[g;b];
end