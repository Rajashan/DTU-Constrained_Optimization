function [x_opt,info] = activeSetQP(H,g,A,b,Aeq,beq)

kmax = 1000;
[varIn,n_cons] = size(A); 
%x_0 = compute_feasiblePhaseOne(A,b,Aeq,beq); %Compute feasible point
x_0 = compute_feasibleBigM(H,g,A,b,Aeq,beq);
%x_0 = compute_feasibleBigMV2(H,g,A,b,Aeq,beq);


if(isempty(x_0))
    error('The algorithm was not able to find a feasible starting point');
end

W_0 = [2];    %Pick random start workingset
workset(1,:) = zeros(1,n_cons);
workset(1,W_0) = 1;

%Define all variables
x(:,1) = x_0;
tol_step = 1e-9;
tol_lambda = 1e-9;
k = 1;
lambda_active = nan(n_cons,1);
while(k < kmax)
    work_idx = find(workset(k,:));
    inactive_set = find(~workset(k,:));
    
    %Solve subproblem to find step-direction
    Aeq_sub = [Aeq, A(:,work_idx)]';     
    beq_sub = zeros(size(Aeq_sub,1),1);
    g_k = H*x(:,k) + g;
    [p(:,k),lam] = EqualityQPSolverActiveSet(H,g_k,Aeq_sub,beq_sub);   
    lambda_active(work_idx,k) = lam;  %Update the lagrange multipliers
    lambda_active(inactive_set,k) = NaN;
    if(norm(p(:,k),'inf') <= tol_step)
        
      %Find the lagrange multipliers
      lambda_active(work_idx,k) = A(:,work_idx)\(H*x(:,k) + g); 
      alpha(k) = NaN;
      if(norm(lambda_active(~isnan(lambda_active(:,k)),k),'inf') >= tol_lambda) %Check for optimum
          x_opt = x(:,k);
          break;
      else
        %Find idx of lowest lagrange multiplier  
        [~,j] = min(lambda_active(:,k)); 
        x(:,k+1) = x(:,k);       
        
        %Remove the constraint with most negative lambda
        workset(k,j) = 0;         
        workset(k+1,:) = workset(k,:); 
        lambda_active(j,k) = NaN;
      end
    else
        inactive_idx = find(~workset(k,:)); %Get idx of inactive constraints
        
        %Calculate 2nd argument
        min_arg2 = (b(inactive_idx) - A(:,inactive_idx)'*x(:,k) ) ./ (A(:,inactive_idx)'*p(:,k)); 
        
        %Finds the idx that fulfills a_i^T*p_k < 0
        valid_min_arg2 = find( (A(:,inactive_idx)'*p(:,k)) < 0);
        alpha(k) = min(1, min(min_arg2(valid_min_arg2)));
        x(:,k+1) = x(:,k) + p(:,k)*alpha(k);
        if(alpha(k) < 1)
          min_arg2(min_arg2 <0) = NaN;
          [~,blockIdx] = min(min_arg2);
          workset(k+1,blockIdx) = 1; 
        else
          workset(k+1,:) = workset(k,:); 
        end
        lambda_active(:,k+1) = lambda_active(:,k); 
        
        
    end    
    k = k +1;
end

info.alpha = alpha;
info.NumIter = k;
info.steps = p;
info.Workingset = workset;
info.lambda = lambda_active;
info.x = x;

end