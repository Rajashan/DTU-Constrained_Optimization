%% Problem 3a - Inequality Constrained SQP
clear variables; close all; clc
options = optimoptions('quadprog','Display','off');

%Starting Guess
%x = [-4 3]';    %The Two Left Boundary Mins
%x = [1 1]';    %Converges to Global Min
%x = [0,0]';     %Converges to Global Min
%x = [3 -2]';    %%Converges to Global Min
%x = [-2 3]';    %The Two Left Boundary Mins
%x = [-0.5 2]';  %Converges to Global Min
y = [0 0]';

%Initialization
maxit = 50;
tol = 10^(-6);
k = 0;
B = 10^2*eye(2);

%Initial Computations
[f df] = ObjFun(x);
[c(1,1) dc(:,1)] = ConFun1(x);
[c(2,1) dc(:,2)] = ConFun2(x);
L = df - dc*y;

%Stats
stat.k = [];
stat.x = x;
stat.y = y;
stat.c = c;
stat.df = df;
stat.dc = dc;
stat.L = L;
stat.B = B;
stat.theta = [];
stat.time = [];

%Main Algorithm
Converged = ( norm(L,'inf') <= tol && ( min(c) >= 0 ) );
while ~Converged && k < maxit
    tic
    k = k + 1
    
    %Solve Inequality Quadratic Program for Step Direction using BFGS
    %Approx by Interior-Point Algorithm (quadprog)
    %%%%%%%%%%%%%%%%%%%%%%%%%
    [p,~,~,~,lambda] = quadprog(B,df,-dc,c,[],[],[],[],[],options);
    y_new = lambda.ineqlin;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %New Step and New Function Evaluations
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    x_new = x + p;
    [f_new df_new] = ObjFun(x_new);
    [c_new(1,1) dc_new(:,1)] = ConFun1(x_new);
    [c_new(2,1) dc_new(:,2)] = ConFun2(x_new);
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
    L = df - dc*y;
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    %Convergence Check
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    Converged = ( norm(L,'inf') <= tol && ( min(c) >= 0 ) );
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Stats Update
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    stat.k = [stat.k k];
    stat.x = [stat.x x];
    stat.y = [stat.y y];
    stat.c = [stat.c c];
    stat.df = [stat.df df];
    stat.dc = [stat.dc dc];
    stat.L = [stat.L L];
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
elseif k >= maxit
    disp('The maximum nubmer of iterations has been reached and no solution was found')
end

%Contour Plot and Iteration Sequence
x = -5:0.05:5;
y = -5:0.05:5;
[X,Y] = meshgrid(x,y);
F = (X.^2+Y-11).^2 + (X + Y.^2 - 7).^2;
v = [0:2:10 10:10:100 100:20:200];
[c,h] = contour(X,Y,F,v,'linewidth',2);
colorbar
yc1 = (x+2).^2;
yc2 = (4*x)/10;
hold on
fill([x(yc1<=5.1)],[yc1(yc1<=5.1)],[0.7 0.7 0.7],'facealpha',0.2)
fill([x x(end) x(1)],[yc2 -5 -5],[0.7 0.7 0.7],'facealpha',0.2)
xlim([-5 5])
ylim([-5 5])

for i = 1:k
    plot(stat.x(1,i),stat.x(2,i),'o','MarkerSize',5,'MarkerEdgeColor','Red','MarkerFaceColor','Red')
    plot(stat.x(1,i:i+1),stat.x(2,i:i+1),'-','Linewidth',1,'Color','Red')
    Numb = sprintf('%d',stat.k(i));
    text(stat.x(1,i)+0.05,stat.x(2,i)+0.05,Numb,'Color','Red','FontSize',10)
end