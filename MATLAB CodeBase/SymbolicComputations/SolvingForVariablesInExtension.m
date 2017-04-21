syms a_1 a_2 a_3 a_4 a_5 a_6 
assume(a_1, 'real')
assume(a_2, 'real')
assume(a_3, 'real')
assume(a_4, 'real')
assume(a_5, 'real')
assume(a_6, 'real')

x = sym('x_',[9,1]);
y = sym('y_',[9,1]);

syms T1 T2 T3
assume(T1, 'real')
assume(T2, 'real')
assume(T3, 'real')

A = toeplitz ( [T1, T2, T3] );
T1 = toeplitz([a_1,a_2, a_3]);
T2 = toeplitz( [a_2 , a_4, a_5]);
T3 = toeplitz([a_3, a_5, a_6]);
A=NumericProdReal;

y = A*x;

syms T1 T2 T3
assume(T1, 'real')
assume(T2, 'real')
assume(T3, 'real')

newA = toeplitz ( [T1, T2, T3] );
T1 = toeplitz([a_3, a_2, a_1]);
T2 = toeplitz( [a_5, a_4, a_2]);
T3 = toeplitz([a_6, a_5, a_3]); 
newA= eval(newA);

y = newA*x;

newA = subs(newA, a_1, NumericProdReal(1,1));
newA = subs(newA, a_2, NumericProdReal(2,1));
newA = subs(newA, a_3, NumericProdReal(3,1) );
newA = subs(newA, a_4, NumericProdReal(5,1) );
newA = subs(newA, a_5, NumericProdReal(6,1) );
newA = subs(newA, a_6, NumericProdReal(9,1) );
newA = eval(newA)


%% A = [1 2 3;4 5 6;7 8 9]; % Example matrix
eig_A = eig(A);
flag = 0;
for i = 1:rank(A)
	if eig_A(i) <= 0 
	flag = 1;
	end
end
if flag == 1
	disp('the matrix A is not positive definite')
	else
	disp('the matrix A is positive definite')
end

eig_A = eig(newA);
flag = 0;
for i = 1:rank(newA)
	if eig_A(i) <= 0 
	flag = 1;
	end
end
if flag == 1
	disp('the matrix newA is not positive definite')
	else
	disp('the matrix newA is positive definite')
end


%%
 
% syms l m n o l_conj m_conj n_conj o_conj
% sym_herm_A = [ 
%          0 0 0 l      1 l_conj 0 0 0 m      1 m_conj; 
%          0 0 0 n      1 n_conj 0 0 0 n_conj 1 n;     
%          0 0 0 m_conj 1 m      0 0 0 l_conj 1 l; 
%          0 0 0 o      1 o_conj 0 0 0 o      1 o_conj; 
%          0 0 0 1      1 1      0 0 0 1      1 1;
%          0 0 0 o_conj 1 o      0 0 0 o_conj 1 o ;
%          0 0 0 m      1 m_conj 0 0 0 l      1 l_conj ; 
%          0 0 0 n_conj 1 n      0 0 0 n      1 n_conj; 
%          0 0 0 l_conj 1 l      0 0 0 m_conj 1 m; 
%          ];
%      
% sym_A = sym_herm_A';
% 
% sym_ProductNew = sym_herm_A* sym_A
% 
% a_1 = sym_ProductNew(1,1);
% a_2 = sym_ProductNew(1,2);
% a_3 = sym_ProductNew(1,3);
% a_4 = sym_ProductNew(1,5);
% a_5 = sym_ProductNew(1,6);
% a_6 = sym_ProductNew(1,9);
% 
% newA = eval(A);
% newA - sym_ProductNew
% subs(newA - sym_ProductNew,l_conj, conj(l))
% subs(newA - sym_ProductNew,m_conj, conj(m))
% subs(newA - sym_ProductNew,o_conj, conj(o))
% subs(newA - sym_ProductNew,n_conj, conj(n))
% % subs(newA - sym_ProductNew,l_conj, conj(l))
% % subs(newA - sym_ProductNew,conj(l_conj), l)
% % subs(newA - sym_ProductNew,m_conj, conj(m))
% % subs(newA - sym_ProductNew,conj(m_conj), m)
% % subs(newA - sym_ProductNew,n_conj, conj(n))
% % subs(newA - sym_ProductNew,conj(n_conj), n)
% % subs(newA - sym_ProductNew,o_conj, conj(o))
% % subs(newA - sym_ProductNew,conj(o_conj), o)

% hermA = A';
% sym_herm_A = hermA;
% sym_herm_A (1,:)= hermA(3,:);
% sym_herm_A (2,:)= hermA(8,:);
% sym_herm_A (3,:)= hermA(1,:);
% sym_herm_A (4,:)= hermA(6,:);
% sym_herm_A (5,:)= hermA(2,:);
% sym_herm_A (6,:)= hermA(4,:);
% sym_herm_A (7,:)= hermA(9,:);
% sym_herm_A (8,:)= hermA(5,:);
% sym_herm_A (9,:)= hermA(7,:);
% 
% 
% sym_ProductNew = sym_herm_A* sym_herm_A';
% sym_ProductNew = subs(sym_ProductNew, alphas, alpha_Values );
% sym_ProductNew = subs(sym_ProductNew,  betas, beta_Values );
% sym_ProductNew = subs(sym_ProductNew, one, 1 );
% sym_ProductNew = eval(sym_ProductNew)
% NumericProdReal

