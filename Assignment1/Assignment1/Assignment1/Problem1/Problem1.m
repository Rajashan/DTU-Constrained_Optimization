%% Exercise 3
clear all; close all;clc
%% Problem 1 - Quadratic Optimization, Sensitivity and Duality

%Function has been created. Let us try call it:
H = [6 2 1; 2 5 2; 1 2 4];
g = -[8 3 3]';
A = [1 0; 0 1; 1 1];
b = [3 0]';

[x,lambda] = EqualityQPSolver(H,g,A,b);

%% Random QPs
n = 3;
m = 2;
[H,g,A,b,sol] = randomQP(n,m);
sol

[x,lambda] = EqualityQPSolver(H,g,A,b);
[x;lambda]

%% 7. Sensitiities

%We have shown on paper that the sensitivities are just inv(K) where K is
%the system matrix

%Lets compute the result for various systems and compare with just the
%sensitivies times the p = [g b]' vector

H = [6 2 1; 2 5 2; 1 2 4];
A = [1 0; 0 1; 1 1];

%We know the sensitivities. Check if two solutions are the same
p = [5 -3 2 8 2]';
g = p(1:3);
b = p(4:5);

y = compSensConvex(H,g,A,b)*p;
[x,lambda] = EqualityQPSolver(H,g,A,b);

%We see that the solutions are indeed the same. So they can be found just
%by the matrix vector product of sensitivities and parameters p. NOTICE
%THIS IS NOT NORMAL BEAVIOUR. USUALLY YOU HAVE TO TAYLOR EXPAND AND USE
%SENSITIVITIES TO APPROXIMATE NEW SOLUTION

%% Sensitivities using Taylor Expansion
%Define the system
H = [6 2 1; 2 5 2; 1 2 4];
A = [1 0; 0 1; 1 1];

%1. Choose p0
p0 = [-8 ; -3 ; -3 ; 3 ; 0];
g0 = p0(1:3);
b0 = p0(4:5);

%2 Choose p0 + epsilon
p = p0 + randi(20,5,1);
g = p(1:3);
b = p(4:5);

%Calculate solution to p0 and p
[x0,lambda0] = EqualityQPSolver(H,g0,A,b0);
[x,lambda] = EqualityQPSolver(H,g,A,b);

%Sensitivities in p0:
dx = compSensConvex(H,g0,A,b0);

%Calculate p solution with taylor expansion
X = [x0;lambda0] + dx'*(p-p0);

%Print comparsion solution
[[x0;lambda0],[x;lambda],X]

%% 8. Solving the Dual Problem

%Specify system matrices
H = [6 2 1; 2 5 2; 1 2 4];
g = -[8 3 3]';
A = [1 0; 0 1; 1 1];
b = [3 0]';

%Construct system
LHS = [H zeros(3,2) -H ; zeros(2,3) zeros(2,2) A' ; -H A zeros(3,3)];
RHS = [zeros(3,1);b;g];

%Solution
sol = LHS \ RHS;

%We see that the solution 1:5 yields the same as the primal problem
sol(1:5)

%The remaining 6:8 are new lagrange multipliers.
sol(6:8)


