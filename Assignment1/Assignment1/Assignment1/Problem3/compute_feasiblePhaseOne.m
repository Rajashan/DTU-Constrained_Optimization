function x_feas = compute_feasiblePhaseOne(A,b,Aeq,beq)

[varIn,n_cons] = size(A);
[~,n_consEq] = size(Aeq);
x_0 = rand(varIn,1);

if(~isempty(Aeq)) %Checks if there is any equality constraints
gamma_eq = -sign(beq - Aeq'*x_0);
Aeq = [Aeq',gamma_eq*eye(n_consEq)];
end

A = [A',eye(n_cons)];
f = [ones(n_cons,1);ones(varIn,1)];   

optimopt = optimoptions('linprog','display','off');
lb = [zeros(1,n_cons), -inf(1,varIn)]; %Only the z's that are bounded to be nonnegative
x_feas = linprog(f,-A,-b,-Aeq,-beq,lb,[],optimopt);
x_feas = x_feas(n_cons+1:end);

end