%% 5 Effective frontier using interior point

H = (1/2).*[2.30, 0.93, 0.62, 0.74, -0.23; 0.93, 1.40, 0.22, 0.56, 0.26;...
    0.62, 0.22, 1.80, 0.78, -0.27; 0.74, 0.56, 0.78, 3.40, -0.56;...
    -0.23, 0.26, -0.27, -0.56, 2.60];

mu = [15.10;12.50;14.70;9.02;17.68];
e1 = ones(5,1);
A = [mu,e1];
g=0;
b = [10;1];

x = ones(size(H,1),1);
y = ones(size(A,2),1);
z = ones(size(x));

C = eye(length(x));
d = zeros(length(x),1);
S = 2.*ones(size(C,1),1);

x_eff = zeros(5,11);
y_eff = zeros(2,11);
z_eff = zeros(5,11);
s_eff = zeros(5,11);
r = [9.02, 9.9,10.7,11.6,12.5,13.3,14.2, 15.0,15.9,16.8,17.6];
var_eff = zeros(1,11);

for i=1:length(r)
    b = [r(i);1];
    [x_eff(:,i),y_eff(:,i),z_eff(:,i),s_eff(:,i),iter]=QPippd10(H,g,C,d,A,b,x,y,z,S);
    var_eff(1,i)=x_eff(:,i)'*H*x_eff(:,i);
end


figure(10)
plot(r,var_eff,'o--')
title('Efficient Frontier')
xlabel('Return (R)')
ylabel('Risk (Variance)')
xlim([min(r),max(r)])
figure(20)
plot(r,x_eff,'-')
title('Portfolio allocation')
xlabel('Return (R)')
ylabel('Proportion')
legend('x_1','x_2','x_3','x_4','x_5')
xlim([min(r),max(r)])

%% 5 Effective frontier using interior point, risk free

H = [2.30, 0.93, 0.62, 0.74, -0.23, 0; 0.93, 1.40, 0.22, 0.56, 0.26, 0;...
    0.62, 0.22, 1.80, 0.78, -0.27, 0; 0.74, 0.56, 0.78, 3.40, -0.56, 0;...
    -0.23, 0.26, -0.27, -0.56, 2.60, 0; 0, 0, 0, 0, 0, 0];
mu = [15.10;12.50;14.70;9.02;17.68;2.0];

e1 = ones(6,1);
A = [mu,e1];
g=0;
%b = [10;1];

x = ones(size(H,1),1);
y = ones(size(A,2),1);
z = ones(size(x));

C = eye(length(x));
d = zeros(length(x),1);
S = 2.*ones(size(C,1),1);

x_eff = zeros(6,11);
y_eff = zeros(2,11);
z_eff = zeros(6,11);
s_eff = zeros(6,11);
r = [2, 3.5,5,6.5,8.0,9.5,11.0, 13.5,15.0,16.5,17.6];
var_eff = zeros(1,11);


for i=1:length(r)
    b = [r(i);1];
    [x_eff(:,i),y_eff(:,i),z_eff(:,i),s_eff(:,i),iter]=QPippd10(H,g,C,d,A,b,x,y,z,S);
    var_eff(1,i)=x_eff(:,i)'*H*x_eff(:,i);
end


figure(10)
plot(r,var_eff,'o--')
title('Efficient Frontier, interior point')
xlabel('Return (R)')
ylabel('Risk (Variance)')
xlim([min(r),max(r)])
figure(20)
plot(r,x_eff,'-')
title('Portfolio allocation, interior point')
xlabel('Return (R)')
ylabel('Proportion')
legend('x_1','x_2','x_3','x_4','x_5','x_6')
xlim([min(r),max(r)])


%% problem from problem 3 MOST RECENT hot start

H = [2 , 0; 0, 2];
C = -[1,-1,-1,1,0;-2,-2,2,0,1];
d = -[-2;-6;-2;0;0];
A = sparse(zeros(size(H)));
b = zeros(length(A(:,1)),1);
g = [-2;-5];

x = ones(2,1);
z = ones(5,1);
y = zeros(2,1);
s = ones(size(z));
e = ones(5,1);
% initial point
rL = H*x + g - C*z;
rC = s - C'*x+d;
rsz = diag(s)*diag(z)*e;

H_bar = H + C*(diag(s)\diag(z))*C';
L = chol(H_bar);

rL_bar = rL-C*(diag(s)\diag(z))*(rC-diag(z)\rsz);

dx_a = L'\(L\-rL_bar);
dz_a = -(diag(s)\diag(z))*C'*dx_a+(diag(s)\diag(z))*(rC-diag(z)\rsz);
ds_a = -diag(z)\(rsz)-(diag(z))\diag(s)*dz_a;

z = max(ones(length(z),1),abs(z+dz_a));
s = max(ones(length(s),1),abs(s+ds_a));


%[x1,y1,z1,s1,iter,seq]=QPippd10(H,g,C,d,A,b,x,y,z,s);
[x_sol,y_sol,z_sol,s_sol,k] = QPippd10(H,g,C,d,A,b,x,y,z,s);
x = quadprog(H,g,C',d,[],[],[],[]);


%% problem from problem 3, cold start
H = [2 , 0; 0, 2];
C = [1,-1,-1,1,0;-2,-2,2,0,1];
d = [-2;-6;-2;0;0];
g = -[2;5];
x = zeros(2,1);
z = 0.1*ones(5,1);
s = 0.1*ones(size(z));
e = ones(5,1);
Z = diag(z);
S = diag(s);


[x2,z2,s2,iter2,seq2]=QPippd11(H,g,C,d,x,z,s);


