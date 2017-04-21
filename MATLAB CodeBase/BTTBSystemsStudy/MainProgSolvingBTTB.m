%% A.8

% The main program for solving BTTB systems T_{mn}u = b. 
clear  
m = input('Input the number of blocks m:');
n = input ('nput the size of each block n:');

disp('choice of preconditioner:')
disp('0: No preconditioner')
disp('1: CB preconditioner')
disp('2: BCCB preconditioner')

pchoice = input ('Input a preconditioner (pchoice): 0, 1, 2:');

fprintf('\n');
disp ('choice of sequences:')
disp('1: Sequence (i) in Table 5.1');
disp('2: Sequence (ii) in Table 5.1');
disp('3: Sequence (iii) in Table 5.2');
disp('4: Sequence (iv) in Table 5.2');

fchoice = input ('choice an example (fchoice): 1^4: ')

ig = zeros(m*n,1);                 % the initial guess
b = ones(m*n, 1);                  % the right-hand side vector
t = fcolrow(m,n, fchoice);         % the first column and first rows
                                   % of each block T_{mn}, see A.9 
tev = gentev(t);                   % the eigenvalues of the BCCB matrix 
                                   % in which T_{mn} is selected , see A.11
tol = 1.e-7;

if pchoice == 2  
    ev = gen12ev(t);
    % call PCG method with the BCCB preconditioner
    u = pcg12(b, ig, tev,ev, tol);
end