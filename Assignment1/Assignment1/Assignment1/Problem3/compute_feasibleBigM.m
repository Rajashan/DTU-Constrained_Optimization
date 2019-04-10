function x_feas = compute_feasibleBigM(H,g,A,b,Aeq,beq)

[varIn,n_cons] = size(A);
[~,n_consEq] = size(Aeq);
M = 2;
eta_tol = 1e-12;
optimopt = optimoptions('quadprog','display','off');



posA = [Aeq',-ones(n_consEq,1)];
posBeq = beq;

negA = -[Aeq',ones(n_consEq,1)];
negBeq = -beq;

inEqA = -[A',ones(n_cons,1)];
inEqb = -b;

A = [posA;negA;inEqA];
b = [posBeq;negBeq;inEqb];

g(end+1) = M;
H(end+1,end+1) = 0;

lb = [-inf(1,varIn), 0];
x_feas = quadprog(H,g,A,b,[],[],lb,[],[],optimopt);

%Solves the problem again with a bigger value of M if eta is larger than tolerance
while(x_feas(end) >= eta_tol)
    M = M^2;                    
    g(end) = M;
    x_feas = quadprog(H,g,A,b,[],[],lb,[],[],optimopt);
end

x_feas = x_feas(1:end-1);

end