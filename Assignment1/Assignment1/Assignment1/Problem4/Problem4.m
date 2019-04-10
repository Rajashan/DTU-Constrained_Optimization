%% 4.3 one solution

H = [2.30, 0.93, 0.62, 0.74, -0.23; 0.93, 1.40, 0.22, 0.56, 0.26;...
    0.62, 0.22, 1.80, 0.78, -0.27; 0.74, 0.56, 0.78, 3.40, -0.56;...
    -0.23, 0.26, -0.27, -0.56, 2.60];

mu = [15.10;12.50;14.70;9.02;17.68];
e = ones(5,1);

A = [mu';e'];
b = [10;1];

C = -eye(5);
d = zeros(size(C,1),1);


[x,var] = quadprog(H,[],C,d,A,b,[],[]);

figure(1)
bar(x)
title('Portfolio allocation for R=10')
xticks([1 2 3 4 5])
xticklabels({'x_1', 'x_2', 'x_3', 'x_4', 'x_5'})
ylabel('Proportion')

%% 4.4 Effective frontier

H = [2.30, 0.93, 0.62, 0.74, -0.23; 0.93, 1.40, 0.22, 0.56, 0.26;...
    0.62, 0.22, 1.80, 0.78, -0.27; 0.74, 0.56, 0.78, 3.40, -0.56;...
    -0.23, 0.26, -0.27, -0.56, 2.60];

mu = [15.10;12.50;14.70;9.02;17.68];
e = ones(5,1);

A = [mu';e'];

C = -eye(5);
d = zeros(size(C,1),1);

x_eff = zeros(5,11);
var_eff = zeros(1,11);
r = [9.02, 9.9,10.7,11.6,12.5,13.3,14.2, 15.0,15.9,16.8,17.6];

for i = 1:length(r)
    b=[r(i),1];
    [x_eff(:,i),var_eff(1,i)]=quadprog(H,[],C,d,A,b,[],[]);
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

%% 4.6 risk free

H_f = [2.30, 0.93, 0.62, 0.74, -0.23, 0; 0.93, 1.40, 0.22, 0.56, 0.26, 0;...
    0.62, 0.22, 1.80, 0.78, -0.27, 0; 0.74, 0.56, 0.78, 3.40, -0.56, 0;...
    -0.23, 0.26, -0.27, -0.56, 2.60, 0; 0, 0, 0, 0, 0, 0];
mu_f = [15.10;12.50;14.70;9.02;17.68;2.0];

e_f = ones(6,1);

A_f = [mu_f';e_f'];


C_f = -eye(6);
d_f = zeros(size(C_f,1),1);

r_f = [2, 3.5,5,6.5,8.0,9.5,11.0, 13.5,15.0,16.5,17.6];

x_efff = zeros(6,length(r_f));
var_efff = zeros(1,length(r_f));
coordinate = zeros(length(H_f(:,1)),2);
for i = 1:length(r_f)
    b_f=[r_f(i),1];
    [x_efff(:,i),var_efff(1,i)]=quadprog(H_f,[],C_f,d_f,A_f,b_f,[],[]);
end

for i=1:length(H_f(:,1))
    coordinate(i,1)=mu_f(i);
    coordinate(i,2)=H_f(i,i);
end

figure(100)
plot(r_f,var_efff,'o--')
title('Efficient Frontier')
xlabel('Return (R)')
ylabel('Risk (Variance)')
xlim([min(r_f),max(r_f)])
hold on
plot(coordinate(:,1),coordinate(:,2),'o')
labels = cellstr(["x_1","x_2","x_3","x_4","x_5","x_6"]);
text(coordinate(:,1),coordinate(:,2),labels)
plot(r,var_eff,'o--')
title('Efficient Frontier, CAL and assets')
xlabel('Return (R)')
ylabel('Risk (Variance)')
xlim([min(r_f),max(r_f)])
legend('Risk free','Efficient frontier')
hold off

figure(200)
plot(r_f,x_efff,'-')
title('Portfolio allocation')
xlabel('Return (R)')
ylabel('Proportion')
legend('x_1','x_2','x_3','x_4','x_5','x_6')
xlim([min(r_f),max(r_f)])


%% 4.4 Effective frontier analytical

H = [2.30, 0.93, 0.62, 0.74, -0.23; 0.93, 1.40, 0.22, 0.56, 0.26;...
    0.62, 0.22, 1.80, 0.78, -0.27; 0.74, 0.56, 0.78, 3.40, -0.56;...
    -0.23, 0.26, -0.27, -0.56, 2.60];

mu = [15.10;12.50;14.70;9.02;17.68];
e = ones(5,1);

A = [mu';e'];

C = -eye(5);
d = zeros(size(C,1),1);


E = inv(H)*A'*inv(A*inv(H)*A');

R =9.02:1:16.02;


y = (A*inv(H)*A')\b;
x = H\A'*y;

figure(20)
%plot(R,0.055*R.^2-1.58*R+12.04)

figure(30)
plot(R,E(:,1).*R+E(:,2))

ylim(0:1)


%analytical
E = inv(H)*A'*inv(A*inv(H)*A');

R =9:1:18;

figure(20)
plot(R,0.055*R.^2-1.58*R+12.04)

figure(30)
plot(R,E(:,1).*R+E(:,2))
ylim([0:1])
