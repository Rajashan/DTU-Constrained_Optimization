function [c,dc,d2c] = ConFun3(x)

c = x(1)^3 + x(2)^3 + 1;

dc = [3*x(1)^2 ; 3*x(2)^2 ; 0 ; 0 ; 0];

d2c = zeros(5,5,1);
d2c(1,1,1) = 6*x(1); d2c(2,2) = 6*x(2);

end