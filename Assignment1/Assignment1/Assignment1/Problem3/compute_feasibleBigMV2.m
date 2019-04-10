function x_feas = compute_feasibleBigMV2(H,g,A,b,Aeq,beq)

M = 10^5;
[varIn,n_cons] = size(A);
[~,n_consEq] = size(Aeq);

Aeq = [Aeq',ones(n_consEq,1),-ones(n_consEq,1)];
A = [A',eye(n_cons)];
lb =  [-inf(1,varIn), zeros(1,2*n_consEq),zeros(1,n_cons)];

%Define the new Hessian and first-order-component vector 
%with extra variables
H(end+1:(end+n_cons+n_consEq),end+1:(end+n_cons+n_consEq)) = 0;
g(end+1:(end+n_cons+n_consEq)) = M;

opt = optimoptions('quadprog','display','off');

x_feas = quadprog(H,g,-A,-b,Aeq,beq,lb,[],[],opt);
x_feas = x_feas(1:end-n_cons);

end