%% Test script for Algorithm 14.3

%Initialize some random variables 

%%% REMEMBER THAT A NEEDS TO HAVE FULL ROW RANK %%%
n = 3;     %Variable
m = 2;     %Constraints


A = 3*randi(4,m,n);   %Generate random A matrix

x = 4*rand(n,1);    %Generate random "solution"
x(m+1:end) = 0;

s = 10*rand(n,1);   %Generate random slack variables
s(1:m) = 0;

b = A*x;            %Calculate equalities from random A and "solution"
lambda = 10*rand(m,1);  %Generate random lagranian multipliers
g = A'*lambda+s;        %Generate obj function

%Check if we get the same as the "solution" 
[x_opt,lambda_opt,s_opt] = PredictorCorrectorV2(g,A,b);