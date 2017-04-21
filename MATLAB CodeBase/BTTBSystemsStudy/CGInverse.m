clc
clear
%% Using cgs with a Matrix Input
A = gallery('wilk',21);
b = sum(A,2);
tol = 1e-12;  maxit = 15; 
M1 = diag([10:-1:1 1 1:10]);
x = cgs(A,b,tol,maxit,M1);

%% Using cgs with a Function Handle
disp(' Now using function handles')
x1 = run_cgs