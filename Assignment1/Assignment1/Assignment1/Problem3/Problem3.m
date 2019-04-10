 %% Problem 3 - Inequality constrained quadratic programming - Example 16.4

%Make a contour plot of the problem

x = linspace(0,5,100);
y = linspace(0,5,100);

[X,Y] = meshgrid(x,y);

Fun = @(x,y) (x-1).^2 + (y - 2.5).^2;
F = Fun(X,Y);
v = [min(min(F)):1:-2, -1. 99:0.2:2, 2.01:max(max(F))];
contour(X,Y,F,v,'linewidth',1.2);
hold on;
yc1 = (x./2 + 1);
yc2 = (-x./2 + 3);
yc3 = (x./2 - 1);
yc3_ny = yc3(yc3>0);

fill([x x(end) x(1)],[yc1 max(y) max(y)], [0.7 0.7 0.7], 'facealpha',0.4);
fill([x x(end) x(1)],[yc2 max(y) max(y)], [0.7 0.7 0.7], 'facealpha',0.4);
fill([x(yc3>0),x(end)] ,[yc3_ny, 0],[0.7 0.7 0.7], 'facealpha',0.4);
xlabel('x_{1}','fontsize',16);
ylabel('x_{2}','fontsize',16);
title('Contour plot of q(x) = (x_{1} - 1)^{2} + (x_{2} - 2.5)^{2}','fontsize',14)
colorbar;
%% Problem 3 - defining the problem

A = [1, -1, -1, 1 ,0; -2 -2 2 0 1];
b = [-2;-6;-2;0;0];
H = 2*eye(2);
g = [-2;-5];
[x_opt,info] = activeSetQP(H,g,A,b,[],[]);

%% Test of feasible point generators
opt = optimoptions('fmincon','algorithm','active-set','StepTolerance',1e-9,'OptimalityTolerance',1e-9,'ConstraintTolerance',1e-10);
fun = @(x) (x(1)-1).^2 + (x(2)-2.5).^2;
x_0 = compute_feasibleBigM(H,g,A,b,[],[]);
x_0v2 = compute_feasibleBigMV2(H,g,A,b,[],[]);

[x_optv1,info] = activeSetQP(H,g,A,b,[],[])
%[x_optv2,info2] = activeSetQP(H,g,A,b,[],[])
[xv1,~,~,outputv1,lambda1] = fmincon(fun,x_0,-A',-b,[],[],[],[],[],opt);
[xv2,~,~,outputv2,lambda2] = fmincon(fun,x_0v2,-A',-b,[],[],[],[],[],opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conceptual QP using Active Set

H = [2 0; 0 2];
g = [-2 -5]';
A = [1 -1 -1 1 0; -2 -2 2 0 1];
b = [-2 -6 -2 0 0]';

d = [0 0 0 0 0]';

%% Iterate 1
xk1 = [2 0]';
gk1 = H*xk1 + g;
W1 = [5,3]';

[p1,lambda1] = EqualityQPSolver(H,gk1,A(:,[3,5]),d([3,5]))
%The solution is p = 0 with negative multipliers. Remove working set {5}
%and try again

xk2 = xk1;
W2 = [5 0]';
%% Iterate 2
[p2,lambda2] = EqualityQPSolver(H,gk1,A(:,[5]),d([5]))

%Solution is x = (-1, 0). Calculate alpha
a1 =  min( 1 , min( ((b([1,4]) - A(:,[1,4]))'*xk2)./((A(:,[1,4]))'*p2) ) )
%We find a = 1

%The next iterate is then xk1 + p2 = (2,0) + (-1,0) = (1,0)
xk3 = (xk2 + p2);
W3 = [5 0]';
%% Iterate 3
gk2 = H*xk3 + g;

[p3,lambda3] = EqualityQPSolver(H,gk2,A(:,[5]),d([5]))

%p3 = 0, and the lagrange multiplier is zero, so we drop W = {5} to the
%empty set

W4 = [0 0]';
xk4 = xk3;

%% Iterate 4
%The empty set has no constraints so A and b must be zero
[p4,lambda4] = EqualityQPSolver(H,gk2,[],[])

%The solution is found to be be p4 = (0 , 2.5)

%Computing alpha (remember only taking the negative A'*p4)
a2 = min( 1 , min( (b([1,2]) - (A(:,[1,2]))'*xk4)./((A(:,[1,2]))'*p4) ) )
%Alpha is 0.6

%Compute the new x value
xk5 = xk4 + p4*a2
% xk3 = xk2 + p4*a2 = (1,0) + 0.6 * (0 , 2.5) = (1 , 1.5)

%This hits constraint 1, so add this to working set
W5 = [1 0]';

%% Iterate 5
gk3 = H*xk5 + g;

[p5,lambda5] = EqualityQPSolver(H,gk3,A(:,1),d(1))
%p = (0.4 , 0.2)

%find alpha
a3 = min( 1 , min( (b(2:3)-(A(:,[2,3]))'*xk5)./(A(:,[2,3]))'*p5 ) )

%Compute new x
xk6 = xk5 + p5 * a3

%Same working set xk4 lies on constraint 1 still
A'*xk6 - b;
W6 = [1 0]';
%% Iterate 6
gk4 = H*xk6 + g;

[p6,lambda6] = EqualityQPSolver(H,gk4,A(:,1),d(1))
%p6 = 0;

%Lambda is positive. We have found the solution to be xk4 = (1.4 , 1.7)





