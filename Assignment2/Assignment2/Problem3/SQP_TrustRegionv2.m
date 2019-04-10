function [x_opt,lambda_opt,info] = SQP_TrustRegionv2(fun,nonLinCons,x)


tol = 1e-6;
eta = 1/4;
gamma = 1/2;
mu = 100;


Trust_reg = 10;

maxIter = 100;
k = 0;
Converged = 0;
lam(:,1) = [1;1];
x(:,1) = x;


while(k < maxIter && ~Converged)
    k = k + 1;
    
    
    
    %%% COMPUTE f, c, GRADIENT OF f %%%
    
    [f_k,g_k,H_k] = fun(x(:,k));
    [c,dc,d2c] = nonLinCons(x(:,k));
    
           
    
    %%% SET UP THE JACOBIAN %%%
    
    A = dc';
        
    %%% CHECK FOR CONVERGENCE %%% 
    if(norm(g_k - A'*lam(:,k),'inf') < tol)
        Converged = 1;
        x_opt = x(:,k);
        lambda_opt = lam(:,k);
        info.x = x;
        info.Iter = k;
        info.lamda = lam;
        info.Steps = p(1:2,:);
    end 
    
    
    %%% SOLVE SUBPROBLEM 18.50 IN N&W %%% 
    A_ineq = [A, eye(2)];
    g_k(end+1:end+2) = mu;
    
    H_k = H_k - lam(1,k)*d2c(:,:,1) - lam(2,k)*d2c(:,:,2);
    xx_L = H_k;
    H_k(end+1:end+2,end+1:end+2) = 0;

    lb = [-Trust_reg, -Trust_reg, 0, 0];
    ub = [Trust_reg,Trust_reg, Inf, Inf];
    
    opt = optimoptions('quadprog','Display','off');
    [p(:,k),~,~,~,lambda] = quadprog(H_k,g_k,-A_ineq,c,[],[],lb,ub,[],opt);
    lam(:,k+1) = lambda.ineqlin;
 
     
    %%% COMPUTE MU %%%
    
    [f_k,df_k,~] = fun(x);
    
    if(p(:,k)' * H_k * p(:,k) > 0)
        sigma = 1;
    else
        sigma = 0;
    end
    
    
    tmp = [df_k;mu;mu]' * p(:,k) + (sigma/2) * (p(:,k)' * H_k * p(:,k) ) / (0.5*norm(c,1));
    if(mu >= tmp)
        mu = mu;
    else
        mu = tmp;
    end    
    
    %Define merit functions from (18.51)
    phi_1 = f_k + mu * sum(max(0,-c)); 
    phi_12 = fun(x(:,k)+p(1:2,k)) + mu * sum(max(0,-nonLinCons(x(:,k)+p(1:2,k)) ) );    
    
    %Define q_mu as in (18.49)
    q_mu = f_k + df_k'*p(1:2,k) + 0.5* (p(1:2,k)' * xx_L * p(1:2,k) ) + mu * sum(max(0,-(c+dc'*p(1:2,k)) ) );
    q_0 = f_k + mu * sum(max(0,-nonLinCons([0;0])) );    
    
    %Calculate rho_k
    rho_k = (phi_1 - phi_12) / (q_0 - q_mu);
    
    
    %%% CHECK IF WE INCREASE REGION OR NOT %%%
    if(rho_k < eta)
        x(:,k+1) = x(:,k) + p(1:2,k);
        Trust_reg = Trust_reg+10;
    else
        x(:,k+1) = x(:,k);
        Trust_reg = gamma * norm(p(1:2,k));
    end
    
    
    
    
end