function [c , dc] = ConFun1(x)
%Function for computing the function value and its gradient at x
x1 = x(1);
x2 = x(2);

c = ((x1+2)^2-x2);

dc = [2*x1+4 ; -1];

end