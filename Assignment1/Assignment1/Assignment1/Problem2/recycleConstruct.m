function [H,g,A,b] = recycleConstruct(n,ubar,d0)

H = 2*eye(n+1);
b = [-d0;zeros(n-1,1)];
g = -ubar*ones(n+1,1);
A = zeros(n+1,n);

for i = 1:n-1
    A(i,i:i+1) = [-1,1];
end

A(n,1) = 1;
A(n:n+1,n) = [-1;-1];

end