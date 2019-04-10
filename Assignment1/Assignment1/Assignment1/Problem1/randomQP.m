function [H,g,A,b,x] = randomQP(n,m)

%n is the number of variables
%m is the number of equalities

% 1. Generate random H, A and solution vector (X,L)
rand_H = randi(20,n,n);
H = rand_H*rand_H';
while(min(eig(H)) <= 0)
    H = H + 0.1*eye(n);
end
A = randi(20,n,m);
X = randi(20,n,1);
L = randi(20,m,1);

%Construct g and b
g = -(H*X-A*L); 
b = (A'*X);
%Concatenate solution
x = [X;L];
end