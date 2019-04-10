function [c,dc,d2c] = evalNonLinConsProb3(x)


c1 = (x(1)+2)^2 - x(2);
c2 = -4*x(1) + 10*x(2);

c = [c1;c2];

dc1_x1 = 2*x(1) + 4;
dc1_x2 = -1;

dc1_grad = [dc1_x1; dc1_x2];

dc2_x1 = -4;
dc2_x2 = 10;

dc2_grad = [dc2_x1; dc2_x2];

dc = [dc1_grad, dc2_grad];

d2c1 = [2 0; 0 0];
d2c2 = [0 0; 0 0];

d2c(:,:,1) = d2c1;
d2c(:,:,2) = d2c2;

end