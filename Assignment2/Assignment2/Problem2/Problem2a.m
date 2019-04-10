%% Problem 2a - Equality Constrained SQP
clear variables; close all; clc

%Starting Guess and 3 Starting Lagrange Multipliers
x = [-1.8 , 1.7 , 1.9 , -0.8 , -0.8]';
y = [ 1 , 1 , 1]';

%Initialization
maxit = 100;
tol = 10^(-14);
k = 0;

%Initial Computations
[f df d2f] = ObjFun(x);
[c(1,1) dc(:,1) d2c(:,:,1)] = ConFun1(x);
[c(2,1) dc(:,2) d2c(:,:,2)] = ConFun2(x);
[c(3,1) dc(:,3) d2c(:,:,3)] = ConFun3(x);

%Stats
stat.k = [];
stat.x = x;
stat.y = y;
stat.df = df;
stat.d2f = d2f;
stat.time = [];

%Main Algorithm
Converged = ( norm(df - dc*y,'inf') <= tol && ( norm(c,'inf') <= tol ) );
while ~Converged && k < maxit
    tic
    k = k + 1;
    %Compute Hessian 
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    G = 0;
    for i = 1:3
        G = G - y(i)*d2c(:,:,i);
    end
    H = d2f + G;
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Solve Quadratic Problem using LDL Factorization
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    KKT = [H , -dc ; -dc' , zeros(3,3)];
    RHS = -[df ; -c];
    
    [L,D,p] = ldl(KKT,'lower','vector');
    d(p) = L' \ ( D \ ( L \ RHS(p) ) ); d = reshape(d,8,1);
    
    dx = d(1:5);
    y = d(6:8);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %New Step and Function Evaluations
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    x = x + dx;
    [f df d2f] = ObjFun(x);
    [c(1,1) dc(:,1) d2c(:,:,1)] = ConFun1(x);
    [c(2,1) dc(:,2) d2c(:,:,2)] = ConFun2(x);
    [c(3,1) dc(:,3) d2c(:,:,3)] = ConFun3(x);
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    %Convergence Check
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    Converged = ( norm(df - dc*y,'inf') <= tol && ( norm(c,'inf') <= tol ) );
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    stat.k = [stat.k k];
    stat.x = [stat.x x];
    stat.y = [stat.y y];
    stat.df = [stat.df df];
    stat.d2f = [stat.d2f d2f];
    stat.time = [stat.time toc];
end

if Converged && k <= maxit
disp('A solution was found and is given by x = ')
disp(x)
disp('The stats are given by')
disp(stat)
end