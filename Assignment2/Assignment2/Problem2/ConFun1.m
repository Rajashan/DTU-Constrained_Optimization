function [c,dc,d2c] = ConFun1(x)

c = x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 + x(5)^2 - 10;

dc = [2*x(1) ; 2*x(2) ; 2*x(3) ; 2*x(4) ; 2*x(5)];

d2c = 2*eye(5);

end