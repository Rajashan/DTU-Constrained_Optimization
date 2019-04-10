function [x,lambda] = recycleNullSpace(n,ubar,d0)

[H,g,A,b] = recycleConstruct(n,ubar,d0);
[Q,Rbar] = qr(A);
m1 = size(Rbar,2);
Q1 = Q(:,1:m1); 
Q2 = Q(:,m1+1:end); 
R = Rbar(1:m1,1:m1);

x_y = (R')\b;

x_z = ((Q2')*H*Q2)\( (-Q2') *(H*Q1*x_y+g));

x  = Q1*x_y + Q2*x_z;

lambda = R\(Q1'*(H*x+g));

end





