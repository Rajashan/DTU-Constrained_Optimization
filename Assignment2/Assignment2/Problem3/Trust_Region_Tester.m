%%
fun = @(x) evalObjFunProb3(x);
nonLin = @(x) evalNonLinConsProb3(x);
x0 = [1.2;4];
[x_opt,lam_opt,info] = SQP_TrustRegionv2(fun,nonLin,x0)

%The SQP_TrustRegion only works for some starting points for some reason
%Some of these are: x = [2;2], [0;1] and [1.2;4]

% Contour plot of problem

x = -5:0.05:5;
y = -5:0.05:5;
[X,Y] = meshgrid(x,y);
F = (X.^2+Y-11).^2 + (X + Y.^2 - 7).^2;
v = [0:2:10 10:10:100 100:20:200];
[c,h]=contour(X,Y,F,v,'linewidth',2);
colorbar
yc1 = (x+2).^2;
yc2 = (4*x)/10;
hold on
fill([x(yc1<=5.1)],[yc1(yc1<=5.1)],[0.7 0.7 0.7],'facealpha',0.2)
fill([x x(end) x(1)],[yc2 -5 -5],[0.7 0.7 0.7],'facealpha',0.2)
plot(x_opt(1,1),x_opt(2,1),'pg','markersize',12,'linewidth',2);
plot(info.x(1,:),info.x(2,:),'ro--','linewidth',1.2);
axis([-5 5 -5 5])