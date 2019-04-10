%% Problem 2 - Test for different solvers

n = 10;
ubar = 0.2;
d0 = 1;
count = 1;
for i = 10:10:1000
   tic; 
   recycleSolverLU(i,ubar,d0);
   t_LU(count) = toc;
   
   tic;
   recycleSolverLDL(i,ubar,d0);
   t_LDL(count) = toc;
   
   tic;
   recycleNullSpace(i,ubar,d0);
   t_NS(count) = toc;
   
   tic;
   recycleRangeSpace(i,ubar,d0);
   t_RS(count) = toc;
   
   n(count) = i;
   count = count + 1;
end

%% Plot the results
plot(n,t_LU,'--','linewidth',1.2)
hold on;
plot(n,t_LDL,'--','linewidth',1.2);
plot(n,t_NS,'--','linewidth',1.2);
plot(n,t_RS,'--','linewidth',1.2);
hold off;
grid on;
xlabel('n','fontsize',14);
ylabel('Time','fontsize',14);
title('CPU-time of each solver as function of n','fontsize',14);
legend({'LU','LDL','Null-Space','Range-Space'},'location','northWest','fontsize',10);

%% Make the matrices sparse now 

sparse_KKT = sparse(recycleKKT(100,ubar,d0));
spy(sparse_KKT);

recycleSparseLUSolver(10,ubar,d0)
recycleSparseLDLSolver(10,ubar,d0)

%% Problem 2 - sparse vs dense CPU-time
count = 1;
for i = 10:10:1000
    tic;
    recycleSparseLUSolver(i,ubar,d0);
    t_LU_sp(count) = toc;
    
    tic;
    recycleSparseLDLSolver(i,ubar,d0);
    t_LDL_sp(count) = toc;
    
    n(count) = i;
    count = count + 1;
end
%% Plot the CPU time for dense vs sparse

plot(n,t_LU,'--','linewidth',1.2)
hold on;
plot(n,t_LDL,'--','linewidth',1.2);
plot(n,t_NS,'--','linewidth',1.2);
plot(n,t_RS,'--','linewidth',1.2);
plot(n,t_LDL_sp,'--','linewidth',1.2);
plot(n,t_LU_sp,'--','linewidth',1.2);
hold off;
grid on;
xlabel('n','fontsize',14);
ylabel('Time','fontsize',14);
title('CPU-time of each solver as function of n','fontsize',14);
legend({'LU','LDL','Null-Space','Range-Space','LDL-sparse','LU-sparse'},'location','northWest','fontsize',10);