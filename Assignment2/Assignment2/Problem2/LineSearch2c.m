function [x_out,alpha] = LineSearch2c(x,y,p,f,df,c)
%Line Search Algorithm for Equality Constrained SQP

%Initial Alpha Guess
alpha = 1;

%Quadratic Approximation Coefficients
b = df'*p - y'*abs(c);
c = f + y'*abs(c);

Converged = false;
while ~Converged
    %Attempted New Step
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    x_new = x + alpha*p;
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Retrieve New Function Evaluations
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    [f_new] = ObjFun(x_new);
    [c_new(1,1)] = ConFun1(x_new);
    [c_new(2,1)] = ConFun2(x_new);
    [c_new(3,1)] = ConFun3(x_new);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Compute q
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    q = f_new + y'*abs(c_new);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Stop or Update Alpha
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    if q <= c + 10^(-1)*b*alpha
        x_out = x_new;
        break
    else
        a = ( q - (c + b*alpha) ) / alpha^2;
        alpha_min = -b/(2*a);
        
        alpha = min( 0.9*alpha , max( alpha_min,0.1*alpha) );
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%% 
end
end