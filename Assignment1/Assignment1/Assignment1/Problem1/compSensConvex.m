function S = compSensConvex(H,g,A,b)
%Computes the Sensitivities for the constrained equility 
%convex quadratic program.

n_constraints = size(A,2);
N = zeros(n_constraints);

K = [H, -A ; -A', N];

%The sensitivities are given by
s = size(K,1);
S = -K\eye(s);

end