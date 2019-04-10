%% Problem 2b - Equality Constrained SQP
clear variables; close all; clc

%Starting Guess and 3 Starting Lagrange Multipliers
x = [-1.8 , 1.7 , 1.9 , -0.8 , -0.8]';
y = [ 1 , 1 , 1]';

%Initialization
maxit = 100;
tol = 10^(-14);
k = 0;
B = eye(5);

%Initial Computations
[f df] = ObjFun(x);
[c(1,1) dc(:,1)] = ConFun1(x);
[c(2,1) dc(:,2)] = ConFun2(x);
[c(3,1) dc(:,3)] = ConFun3(x);

%Stats
stat.k = [];
stat.x = x;
stat.y = y;
stat.df = df;
stat.B = B;
stat.theta = [];
stat.time = [];

%Main Algorithm
Converged = ( norm(df - dc*y,'inf') <= tol && ( norm(c,'inf') <= tol ) );
while ~Converged && k < maxit
    tic
    k = k + 1;
    
    %Solve BFGS Problem using LDL Factorization
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    KKT = [B , -dc ; -dc' , zeros(3,3)];
    RHS = -[df ; -c];
    
    [L,D,s] = ldl(KKT,'lower','vector');
    d(s) = L' \ ( D \ ( L \ RHS(s) ) ); d = reshape(d,8,1);
    
    p = d(1:5);
    y_new = d(6:8);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %New Step and New Function Evaluations
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    x_new = x + p;
    [f_new df_new] = ObjFun(x_new);
    [c_new(1,1) dc_new(:,1)] = ConFun1(x_new);
    [c_new(2,1) dc_new(:,2)] = ConFun2(x_new);
    [c_new(3,1) dc_new(:,3)] = ConFun3(x_new);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %BFGS Update Step
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    q = ( df_new - dc_new*y_new ) - ( df - dc*y_new ) ;
    
    if p'*q >= 0.2 * p'*(B*p)
        theta = 1;
    else
        theta =  ( 0.8*p'*(B*p) ) / ( p' * (B*p) - p'*q );
    end
    
    r = theta * q + (1-theta)*(B*p);
    
    B = B + (r * r') ./ ( p'*r ) - (B*p)*(B*p)' / ( p'*(B*p) );
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Rename Variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    f = f_new; df = df_new;
    c = c_new; dc = dc_new;
    y = y_new;
    x = x_new;
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Convergence Check
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    Converged = ( norm(df - dc*y,'inf') <= tol && ( norm(c,'inf') <= tol ) );
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Stats Update
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    stat.k = [stat.k k];
    stat.x = [stat.x x];
    stat.y = [stat.y y];
    stat.df = [stat.df df];
    stat.B = [stat.B B];
    stat.theta = [stat.theta theta];
    stat.time = [stat.time toc];
    %%%%%%%%%%%%%%%%%%%%%%%%%%
end

if Converged && k <= maxit
disp('A solution was found and is given by x = ')
disp(x)
disp('The stats are given by')
disp(stat)
end