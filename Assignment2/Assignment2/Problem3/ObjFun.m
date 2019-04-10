function [f , df] = ObjFun(x)
%Function for computing the function value and its gradient at x
x1 = x(1);
x2 = x(2);

f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2;

df = [4*x1^3+2*x2^2+4*x1*x2-42*x1-14 ; 4*x2^3+2*x1^2+4*x1*x2-26*x2-22];

end