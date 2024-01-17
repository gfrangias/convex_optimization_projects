clc, close all, clear all

%%%%%%%%%%%%%%%%%B.i.(a)%%%%%%%%%%%%%%%%%%%
n=2;
A = rand(n,n);

%clc, fprintf('\n Random matrix A is:\n')
%disp(A), fprintf('\n Press a key to continue...'), pause

[U,S,V] = svd(A);
U_1=U*U';
U_2=U'*U;

%clc, fprintf('\n U*U'' is:\n')
%disp(U_1), fprintf('\n Press a key to continue...'), pause
%clc, fprintf('\n U''*U is:\n')
%disp(U_2), fprintf('\n Press a key to continue...'), pause

%%%%%%%%%%%%%%%%%B.i.(b)%%%%%%%%%%%%%%%%%%%

K = 1000;
l_min = abs(randn(1,1));
l_max = K*l_min;

z = l_min + (l_max-l_min)*rand(n-2,1);
eig_P = [l_min;l_max;z];
Lambda = diag(eig_P);

%%%%%%%%%%%%%%%%%%B.ii.%%%%%%%%%%%%%%%%%%%%
q = rand(n,1);
P = U * Lambda * U';

%%%%%%%%%%%%%%%%%%B.iii.%%%%%%%%%%%%%%%%%%%%
x_opt = linsolve(P,-q);
f_opt = 0.5*x_opt'*P*x_opt+q'*x_opt;

fprintf('\n Matrix P is:\n')
disp(P),
fprintf('\n The optimal x is:\n')
disp(x_opt),
fprintf('\n The optimal value of the function using closed-form solution is: %4.2f\n',f_opt),
fprintf('\n Press a key to continue...\n'),
pause

%%%%%%%%%%%%%%%%%%B.iv.%%%%%%%%%%%%%%%%%%%%
cvx_begin
    variable x(n)
    minimize( 0.5*x'*P*x+q'*x )
cvx_end
fprintf('\n Press a key to continue...\n\n'),
pause

%%%%%%%%%%%%%%%%%%B.v.%%%%%%%%%%%%%%%%%%%%
x0=rand(n,1); 
s=1;        % always 1
alpha=0.5;  % from 0 to 0.5
beta=0.1;   % from 0 to 1
epsilon=0.1;
f = @(x) 0.5*x'*P*x+q'*x;
g = @(x) P*x+q;

fprintf('\n Exact Line Search ...\n\n');
[x_e,fun_val_e,iter_e,x_hist_e,fun_val_hist_e] = ...
    gradient_method_exact(f,g,x0,P,epsilon);
fprintf('\n Exact Line Search ^\n\n');

fprintf('\n Press a key to continue...\n\n'),
pause

fprintf('\n Backtracking Line Search ...\n\n');
[x_b,fun_val_b,iter_b,x_hist_b,fun_val_hist_b] = ...
    gradient_method_backtracking(f,g,x0,s,alpha,beta,epsilon);
fprintf('\n Backtracking Line Search ^\n\n');

%%%%%%%%%%%%%%%%%%B.vi.%%%%%%%%%%%%%%%%%%%%
if(n==2)
    max_dist = max(abs(x_opt(1,:)-x0(1)),abs(x_opt(2,:)-x0(2)));
    x1 = (x_opt(1,:)-max_dist-1:0.01:x_opt(1,:)+max_dist)+1;
    x2 = (x_opt(2,:)-max_dist-1:0.01:x_opt(2,:)+max_dist)+1;
    for ii=1:length(x1)
        for jj=1:length(x2)
            x3 = [x1(ii);x2(jj)];
            f2(ii,jj) = 0.5*x3'*P*x3+q'*x3;
        end
    end
    
    figure(1)
    subplot(1,2,1);
    hold on;
    contour(x1, x2, f2'), xlabel('x1'), ylabel('x2'), title('Level sets of the function')
    plot(x_hist_e(1,:),x_hist_e(2,:),"b.-")
    subplot(1,2,2);
    hold on;
    contour(x1, x2, f2'), xlabel('x1'), ylabel('x2'), title('Level sets of the function')
    plot(x_hist_b(1,:),x_hist_b(2,:),"r.-")

    quantity_e = log(fun_val_hist_e-f_opt);
    k_e = (0:iter_e);
    quantity_b = log(fun_val_hist_b-f_opt);
    k_b = (0:iter_b);
    
    %%%%%%%%%%%%%%%%%%B.vii.%%%%%%%%%%%%%%%%%%%%
    figure(2)
    subplot(1,2,1);
    plot(k_e,quantity_e,"b.-"), xlabel('k'), ylabel('log(f(x_k)-p_\ast)')
    subplot(1,2,2);
    plot(k_b,quantity_b,"r.-"), xlabel('k'), ylabel('log(f(x_k)-p_\ast)')
    clc, fprintf('\nI plot the level sets of the function. Press a key to continue...\n'), pause

    %%%%%%%%%%%%%%%%%B.viii.%%%%%%%%%%%%%%%%%%%

    k_max = K*log((fun_val_hist_e(1)-f_opt)/epsilon);
    
    clc, fprintf('\n Maximum Νumber of Ιterations: %0.f\n', k_max),
    fprintf('\n Number of iterations:\n'),
    fprintf(' Exact Line Search: %3d\n',iter_e),
    fprintf(' Backtracking Line Search: %3d\n',iter_b)
    

end