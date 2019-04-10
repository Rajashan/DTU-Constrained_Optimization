function [f,df,d2f] = ObjFun(x)
%Computation of the function its derivative and its hessian

%Variables
x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4); x5 = x(5);

%Useful expressions
e = exp(x1*x2*x3*x4*x5);
k = (x1^3 + x2^3 + 1);

%Function value
f = e - 0.5*k^2;

%Derivative
df = [ x(2)*x(3)*x(4)*x(5)*e - 3*k*x(1)^2 ; x(1)*x(3)*x(4)*x(5)*e - 3*k*x(2)^2 ; x(1)*x(2)*x(4)*x(5)*e ; x(1)*x(2)*x(3)*x(5)*e ; x(1)*x(2)*x(3)*x(4)*e ];

%Hessian

H1 = [ x2^2*x3^2*x4^2*x5^2*e - 9*x1^4 - 6*k*x1 , x3*x4*x5*e + x1*x2*x3^2*x4^2*x5^2*e - 9*x2^2*x1^2 , x2*x4*x5*e + x1*x2^2*x3*x4^2*x5^2*e , x2*x3*x5*e + x1*x2^2*x3^2*x4*x5*e , x2*x3*x4*e + x1*x2^2*x3^2*x4^2*x5*e ];
H2 = [ x3*x4*x5*e + x1*x2*x3^2*x4^2*x5^2*e - 9*x1^2*x2^2 , (x1*x2*x3*x4)^2*e - 9*x2^4 - 6*k*x2 , x1*x4*x5*e + x1^2*x2*x3*x4^2*x5^2*e , x1*x3*x5*e + x1^2*x2*x3^2*x4*x5^2*e , x1*x3*x4*e + x1^2*x2*x3^2*x4^2*x5*e ];
H3 = [ x2*x4*x5*e + x1*x2^2*x3*x4^2*x5^2*e , x1*x4*x5*e + x1^2*x2*x3*x4^2*x5^2*e , x1^2*x2^2*x4^2*x5^2*e , x1*x2*x5*e + x1^2*x2^2*x3*x4^2*x5^2*e , x1*x2*x4*e + x1^2*x2^2*x3*x4^2*x5*e ];
H4 = [ x2*x3*x5*e + x1*x2^2*x3^2*x4*x5^2*e , x1*x3*x5*e + x1^2*x2*x3^2*x4*x5^2*e , x1*x2*x5*e + x1^2*x2^2*x3*x4*x5^2*e , (x1*x2*x3*x5)^2*e , x1*x2*x3*e + x1^2*x2^2*x3^2*x4*x5*e ];
H5 = [ x2*x3*x4*e + x1*x2^2*x3^2*x4^2*x5*e , x1*x3*x4*e + x1^2*x2*x3^2*x4^2*x5*e , x1*x2*x4*e + x1^2*x2^2*x3*x4^2*x5*e , x1*x2*x3*e + x1^2*x2^2*x3^2*x4*x5*e , (x1*x2*x3*x4)^2*e ];
d2f = [H1 ; H2 ; H3 ; H4 ; H5];


end