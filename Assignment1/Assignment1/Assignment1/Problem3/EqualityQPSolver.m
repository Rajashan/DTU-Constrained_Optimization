function [x,lambda] = EqualityQPSolver(H,g,A,b)
%Solves an equality constrained quadratic programme. Specifically written
%for Exercise3

%We require the number of constraints to generate a zero matrix
n_constraints = size(A,2);
N = zeros(n_constraints);

%No. of Dimensions
dim = size(H,1);

%The system is given by
LHS = [H, -A ; -A', N];
RHS = -[g ; b];

%Solution. Should be factorized using LDL since LHS is symmetric
% and indefinite in the general case. Especially useful for many iterates
% of different right-hand sides.
[L,D,p] = ldl(LHS,'lower','vector');
sol(p) = L' \ ( D \ ( L \ RHS(p) ) );

%Exctracting x and lambda from solution vector
x = sol(1:dim)';
lambda = sol(dim+1:end)';
end

