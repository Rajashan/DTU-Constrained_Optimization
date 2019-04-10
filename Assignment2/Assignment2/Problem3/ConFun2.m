function [c , dc] = ConFun2(x)
%Function for computing the function value and its gradient at x
x1 = x(1);
x2 = x(2);

c = -4*x1+10*x2;

dc = [-4 ; 10];

end