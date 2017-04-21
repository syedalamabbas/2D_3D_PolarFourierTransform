clc

% n = 7;    % Must be odd
% [I,J] = ndgrid(1:n);
% A = mod(I+J+(n-3)/2,n);
% B = mod(I+2*J-2,n);
% M = n*A + B + 1
% magic(n) 
% disp('MATLAB algo')
% disp(magic(n))
% syms m
% det(M -m*eye(n) )
% trace(M) 
% eig(M)
% 
% eig((n-1)/2*A + B + 1)
% eig((n+1)/2*A + B + 1)


%% All symbolic
% a = sym('a', [ n 1 ])
% 
% I = a*ones(1,n)
% J = ones(n,1)*a'
% A = mod(I+J+(n-3)/2,n)
% B = mod(I+2*J-2,n)
% M = n*A + B + 1
% tr_M = trace(M) 
% 
% tr_M = subs(tr_M, a, [1:n]')

% 
% n = 132; 
% k = (-n/2:n/2-1);
% l = (-n/2:n/2-1);
% u =(-n/2:n/2-1);
% thier_Toeplitz = zeros(n, n);
% 
%  for x = 1:n          % index for k 
%      for y = 1:n       % index for l 
%          for z=1:n      % index for u 
%              thier_Toeplitz(x,y) = thier_Toeplitz(x,y) + exp(4*pi*1i*u(z)*(k(x)-l(y))/ (2*n+1));
%          end
%      end
%  end
%  thier_Toeplitz
%  cond(thier_Toeplitz)
 

syms a_1 a_2 a_3 a_4 a_5 a_6
assume(a_1, 'real')
assume(a_2, 'real')
assume(a_3, 'real')
assume(a_4, 'real')
assume(a_5, 'real')
assume(a_6, 'real')
N = 2;

%% Constructing a full BTTB Matrix
BTTB_ColumnsIn2DShape = [a_1 a_2 a_3; a_2 a_4 a_5; a_3 a_5 a_6];
Indexes = 1:(N+1);
Toeplitz_Indexes = toeplitz(Indexes);   % Block Level

BTTB_Matrix = []; 
for l = 1:(N+1)
    ColumnIndexes = Toeplitz_Indexes(:,l);
    ToeplitzBlocks = [];
    for k = 1:(N+1)
        CurrentBlockIndex = ColumnIndexes(k);                      % This is the block Index
        ToeplitzBlock = toeplitz(BTTB_ColumnsIn2DShape(:,CurrentBlockIndex));
        ToeplitzBlocks = [ ToeplitzBlocks ; ToeplitzBlock];    % Append Column-wise T = [T1 ; T2; T3]
    end
    BTTB_Matrix  = [BTTB_Matrix, ToeplitzBlocks];              % Append Row-wise
end

BTTB_Matrix


%% Constructing a full BCCB matrix

