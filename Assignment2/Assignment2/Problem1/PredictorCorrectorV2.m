function [x,lambda,s] = PredictorCorrectorV2(g,A,b)


[m,n] = size(A);

eta = 0.995;
k = 0;
kmax = 100;
tolL = 1.0e-9;
tolA = 1.0e-9;
tols = 1.0e-9;

%Calculate x_0, lambda_0 and s_0 as described on page 410 N&W
len = length(A*A');
invA = (A*A')\eye(len);
x_tilde = A'*invA*b;
lambda_tilde = invA*A*g;
s_tilde = g - A'*lambda_tilde;

delta_x = max(-3/2,min(min(x_tilde,0)));
delta_s = max(-3/2,min(min(s_tilde,0)));

x_tilde = x_tilde + delta_x*ones(length(x_tilde),1);
s_tilde = s_tilde + delta_s*ones(length(s_tilde),1);

delta_x2 = 0.5 * (x_tilde'*s_tilde)/(ones(1,length(s_tilde))*s_tilde);
delta_s2 = 0.5 * (x_tilde'*s_tilde)/(ones(1,length(x_tilde))*x_tilde);

x = x_tilde + delta_x2*ones(length(x_tilde),1);
s = s_tilde + delta_s2*ones(length(s_tilde),1);
lambda = lambda_tilde;


%Check residuals 
r_c = -(A'*lambda + s - g); %Dual constraints
r_b = A*x - b;              %Primal constraints
XSe = x.*s;                 %Complementarity conditions
mu  = sum(XSe)/n;           %Duality measure

% Check if initial guess is optimal Converged
Converged = (norm(r_c,inf) <= tolL) && ...
            (norm(r_b,inf) <= tolA) && ...
            (abs(mu) <= tols);



while(k < kmax && ~Converged)
   
   k = k + 1; 
    
   %%%%   AFFINE STEP %%%%%%% 
   %Define left hand side for computing affine direction as in Lecture 8
   lhs_affine = [zeros(n), A', eye(n); A, zeros(m,m), zeros(m,n);...
                 diag(s), zeros(n,m), diag(x)];
   
   %%% FACTORIZATION %%%%
   [L,U] = lu(lhs_affine);
             
   %Define right hand side for computing affine direction as in Lecture 8          
   rhs_affine  = [-r_c; -r_b; -XSe];
   
   %Solve system for affine direction
   affine_direction = U\(L\rhs_affine);
   
   %Extract the different values for each parameter
   x_affine(:,1)      = affine_direction(1:n);
   lambda_affine(:,1) = affine_direction(n+1:end-n);
   s_affine(:,1)      = affine_direction(end-n+1:end);
   
   %%%%% STEP LENGTH %%%%%
   idx = x_affine < 0;
   alpha_primal_aff   = min([1; -x(idx,1)./x_affine(idx,1)]);
   
   idx = s_affine < 0;
   alpha_dual_aff     = min([1; -s(idx,1)./s_affine(idx,1)]);
   
   %%%% CENTERING PARAMETER %%%%%
   xAff_step = x + alpha_primal_aff * x_affine;
   sAff_step = s + alpha_dual_aff * s_affine;
   
   mu_Aff = (xAff_step' * sAff_step)/n;
   
   sigma = (mu_Aff/mu)^3;
   
   
   %%%% SEARCH DIRECTION (14.35) in N&W %%%%
   
   rhs = [-r_c; -r_b; (-XSe -(x_affine .* s_affine)) + sigma*mu];
   
   search_Direction = U\(L\rhs);
   dx(:,1)       = search_Direction(1:n);
   dLambda(:,1)  = search_Direction(n+1 : end-n);
   ds(:,1)       = search_Direction(end-n+1 : end);
   
   %%%%% TAKE THE FINAL STEP %%%%%%
   
   %%% Find step length as in (14.38) %%%
   idx = dx < 0;   
   alpha_pri = min([1; eta * (-x(idx,1)./dx(idx,1))]);
    
   idx = ds < 0;
   alpha_dual = min([1; eta * (-s(idx,1)./ds(idx,1))]);
   
   %%% Compute the step %%% 
   
   x = x + alpha_pri * dx;
   lambda = lambda + alpha_dual * dLambda;
   s = s + alpha_dual * ds;
   
   %%%% UPDATE RESIDUALS  %%%%%
   r_c = A'*lambda + s - g;
   r_b = A*x - b;
   XSe = x.*s;
   mu  = sum(XSe)/n;          

% Converged
Converged = (norm(r_c,inf) <= tolL) && ...
            (norm(r_b,inf) <= tolA) && ...
            (abs(mu) <= tols);

   
end

end