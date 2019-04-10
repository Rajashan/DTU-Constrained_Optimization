function [x,lambda] = EqualityQPSolverActiveSet(H,g,A,b)

rhs = -[g;b];
lhs = [H, -A'; -A zeros(size(A',2))];


[L,D,p] = ldl(lhs,'lower','vector');  %Use LDL factorizatin because Matrix is symmetric and indefinite

sol(p) = L'\( D \(L\rhs(p)));



x = sol(1:length(H));
lambda = sol(length(H)+1:end);

x = x';
lambda = lambda';

end