function [f,df,d2f] = evalObjFunProb3(x)

%Evaluate objective function
f = (x(1)^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2;

%Evaluate gradient of objective function
df1 = 4*x(1)^3 + 4*x(1)*x(2) + 2*x(2)^2 - 42*x(1) - 14;
df2 = 4*x(2)^3 + 2*x(1)^2 + 4*x(1)*x(2) - 26*x(2) - 22;

df = [df1;df2];


%Evaluate hessian of objective function
d2f = [12*x(1)^2 + 4*x(2) - 42, 4*x(1)+4*x(2); 4*x(1)+4*x(2),...
       12*x(2)^2 + 4*x(1) - 26]; 

end