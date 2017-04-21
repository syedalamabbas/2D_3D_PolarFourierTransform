clc
clear
close all

N = 7;   % always odd

x = sym('x_', [1,N]);
syms a CenteredFrFT
assume (a, 'real')

CenteredFrFT =  sym ('null',[N,N]);

middleValue = (N-1)/2;
initialIndex = middleValue+1;

unaliasedArray = -middleValue:middleValue;

for g = 1: N
    CenteredFrFT(g, :) = exp (-1i* 2*pi* unaliasedArray(g) * a * unaliasedArray' /N);
end

disp('This is the FrFT matrix')
CenteredFrFT

disp('Running LU decomposition')
[L, U] = lu(CenteredFrFT);

L = simplify(L, 'IgnoreAnalyticConstraints', true);
simplify((L), 'IgnoreAnalyticConstraints', true)
U = simplify(U,  'IgnoreAnalyticConstraints', true);
simplify((U),  'IgnoreAnalyticConstraints', true)

%  disp('Running svd decomposition')
%  simplify( eig(CenteredFrFT) )