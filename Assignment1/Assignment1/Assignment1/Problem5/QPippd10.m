function [x_sol,y_sol,z_sol,s_sol,k] = QPippd10(H,g,C,d,A,b,x,y,z,s)
% LPIPPD   Primal-Dual Interior-Point QP Solver
%
%          min  x'Hx+g'x
%           x
%          s.t. A x  = b      
%                 C x >= d      
%
% Syntax: [x_s,y_s,z_s,s_s,k] = QPIPPD(H,g,C,d,A,b,x,y,z,S)
%

eta = 0.995;

%residuals. Dimensions are not all used,
% useful for checking

[mA,nA]=size(A);
[mC,nC]=size(C);
e = ones(nC,1);
Z = diag(z);
S = diag(s);


rL = H*x+g-A*y-C*z;
rA = -A'*x+b;
rC = -C'*x+s+d;
rSZ = S*Z*e;
mu = z'*s/nC;

% iteration stopping criteria
k = 0;
maxit = 100;
tolL = 1.0e-9;
tolA = 1.0e-9;
tolC = 1.0e-9;
tolmu = 1.0e-9;



while (k<=maxit && norm(rL)>=tolL && norm(rA)>=tolA && norm(rC)>=tolC ...
        && abs(mu)>=tolmu)
    % factorization of the lhs using LDL  
    H_bar = H + C*(S\Z)*C';
    
    lhs = [H_bar, -A;-A',zeros(size(A,2))];
    [L,D,P] = ldl(lhs,'lower','vector');
 
    % affine direction
    
    rL_bar = rL-C*(S\Z)*(rC-Z\rSZ);
    rhs = -[rL_bar;rA];
    dxy_a(P,:) = L'\(D\(L\(rhs(P,:))));
    
    dx_a = dxy_a(1:length(x));
    dy_a = dxy_a(length(x)+1:length(x)+length(y));
    
    dz_a = -(S\Z)*C'*dx_a+(S\Z)*(rC-Z\rSZ);
    ds_a = -(Z\rSZ)-(Z\(S*dz_a));
    
    % compute largest alpha so we preserve complementarity
    
    alpha_a = 1;
    idx_z = find(dz_a<0);
    if (isempty(idx_z)==0)
        alpha_a = min(alpha_a,min(-z(idx_z)./dz_a(idx_z)));
    end
    idx_s = find(ds_a<0);
    if (isempty(idx_s)==0)
        alpha_a = min(alpha_a,min(-s(idx_s)./ds_a(idx_s)));
    end
    
    
    % affine duality gap
    
    mu_a = ((z+alpha_a*dz_a)'*(s+alpha_a*ds_a))/nC;
    
    % centering parameter (conventional choice)
    
    sigma = (mu_a/mu)^3;
    
    % corrector, using same factorization of lhs
   
    rSZ_bar = rSZ + diag(ds_a)*diag(dz_a)*e - sigma*mu*e;
    rL_bar = rL-C*(S\Z)*(rC-Z\rSZ_bar);
    rhs = -[rL_bar;rA];
    dxy(P,:) = (L'\(D\(L\(rhs(P,:)))));
    dx = dxy(1:length(x));
    dy = dxy(length(x)+1:length(x)+length(y));
    
    dz = -(S\Z)*C'*dx+(S\Z)*(rC-Z\rSZ_bar);
    ds = -Z\rSZ_bar-Z\S*dz;
    
   
    % alpha in the corrector step
    
    alpha = 1;
    idx_z = find(dz<0);
    if (isempty(idx_z)==0)
        alpha = min(alpha,min(-z(idx_z)./dz(idx_z)));
    end
    idx_s = find(ds<0);
    if (isempty(idx_s)==0)
        alpha = min(alpha,min(-s(idx_s)./ds(idx_s)));
    end
    
    % update solutions
    
    x = x + eta*alpha*dx;
    y = y + eta*alpha*dy;
    z = z + eta*alpha*dz;
    s = s + eta*alpha*ds;
    
    Z = diag(z);
    S = diag(s);
    k = k + 1;
    
    rL = H*x + g - A*y- C*z;
    rA = -A'*x+b;
    rC = -C'*x + s + d;
    rSZ = S*Z*e;
    mu =z'*s/nC;
    
end

x_sol = x;
y_sol = y;
z_sol = z;
s_sol = s;
    
    