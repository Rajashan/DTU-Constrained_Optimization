function [c,dc,d2c] = ConFun2(x)

c = x(2)*x(3) - 5*x(4)*x(5);

dc = [0 ; x(3) ; x(2) ; -5*x(5) ; -5*x(4)];

d2c = zeros(5);
d2c(3,2) = 1; d2c(2,3) = 1; d2c(5,4) = -5; d2c(4,5) = -5;

end