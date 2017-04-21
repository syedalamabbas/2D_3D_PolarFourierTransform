clc
clear
close all


%% Making a symmetric Toeplitz into Circulant
N = 2;  % Always even
x = sym('x_', [  (N+1)^2, 1]);
assume (x,'real')


padded_x = padarray(zeros (N+1,N+1), [ N-1,N-1],'post');
sym_padded_x = sym(padded_x);  % converting to symbolic Array
sym_padded_x (1:N+1, 1:N+1) = reshape(x, [N+1,N+1]);
sym_padded_x = reshape(sym_padded_x , [(2*N)^2,1]);

syms x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8 x_9 x_10 x_11 x_12 x_13 x_14 x_15
assume(x_1, 'real')
assume(x_2, 'real')
assume(x_3, 'real')
assume(x_4, 'real')
assume(x_5, 'real')
assume(x_6, 'real')
assume(x_7, 'real')
assume(x_8, 'real')
assume(x_9, 'real')
assume(x_10, 'real')
assume(x_11, 'real')
assume(x_12, 'real')
assume(x_13, 'real')
assume(x_14, 'real')
assume(x_15, 'real')


% x = [x_1; x_2; x_3; x_4; x_5; x_6; x_7; x_8; x_9]  
% 
% 
% y = sym('y_', [  (N+1)^2, 1]);
% assume (y,'real')
% y
% 
% syms y_1 y_2 y_3 y_4 y_5 y_6 y_7 y_8 y_9

% padded_x =  [x(1:3)',0,x(4:6)',0, x(7:9)',0, 0,0,0,0 ]';                 % hat_x

% padded_x =  [x(1:3)',x(2),x(4:6)',x(5), x(7:9)',x(8), 0,0,0,0 ]';                 % hat_x

% padded_y =  [y(1:3)',y(2),y(4:6)',0, y(7:9)',0, 0,0,0,0 ]';                 % hat_y
% padded_x
%% First making just BTTB matrix
% Making a BCCB matrix out of a column of BTTB
% BTTB_ColumnsIn2DShape = sym('a_', [N+1, (N+1)]);
% assume(BTTB_ColumnsIn2DShape, 'real');
syms a_1 a_2 a_3 a_4 a_5 a_6
assume(a_1, 'real')
assume(a_2, 'real')
assume(a_3, 'real')
assume(a_4, 'real')
assume(a_5, 'real')
assume(a_6, 'real')


syms a_11 a_21 a_31 a_12 a_22 a_32 a_13 a_23 a_33
assume(a_11, 'real')
assume(a_21, 'real')
assume(a_31, 'real')
assume(a_12, 'real')
assume(a_22, 'real')
assume(a_32, 'real')
assume(a_13, 'real')
assume(a_23, 'real')
assume(a_33, 'real')

if (N == 2)
    
%     BTTB_ColumnsIn2DShape = [a_1 a_2 a_3; a_2 a_4 a_5; a_3 a_5 a_6];       % N = 2
     BTTB_ColumnsIn2DShape = [a_11 a_21 a_31; a_12 a_22 a_32; a_13 a_23 a_33]';       % N = 2
end

if(N == 4)
    syms a_7 a_8 a_9 a_10 a_11 a_12 a_13 a_14 a_15
    assume(a_7, 'real')
    assume(a_8, 'real')
    assume(a_9, 'real')
    assume(a_10, 'real')
    assume(a_11, 'real')
    assume(a_12, 'real')
    assume(a_13, 'real')
    assume(a_14, 'real')
    assume(a_15, 'real')
    
    BTTB_ColumnsIn2DShape =  [a_1 a_2 a_3 a_4 a_5; a_2 a_6 a_7 a_8 a_9; a_3 a_7 a_10 a_11 a_12; a_4 a_8 a_11 a_13 a_14; a_5 a_9 a_12 a_14 a_15];       % N = 4
end


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



% BTTB_Matrix * x

% [L,U] = lu(BTTB_Matrix )
% % y = BTTB_Matrix * x;
% % disp(y) 
% 
% simplify(y(1)- y(3))
% simplify(y(1)- y(7))
%  
% [X, R] = linsolve(BTTB_Matrix, y)
% 
% 
% hat_y4 = y_2 + (a_3-a_1)*X(2)+ (a_5-a_2)*X(5)+ (a_6-a_3)*X(8);
% 
% hat_y4 = simplify(hat_y4 );
% s = subs (hat_y4,[ y_2, y_3, y_4, y_5, y_6, y_7, y_8, y_9], [0,0,0,0,0,0,0,0] )
% s1 =  subs (s, y_1,1)

%% Making just the BCCB matrix
BCCB_Column = [];                                 % Making each block circulant
for k= 1:N+1
    bttb_col =  BTTB_ColumnsIn2DShape(:,k);
    flippedColumn = flipud(bttb_col);
    BCCB_Column = [BCCB_Column; bttb_col; flippedColumn(2:N) ];
end

flippedIndexes = fliplr(2:N);                   % Making full matrix block circulant
for k= 1:N-1
    bttb_col = BTTB_ColumnsIn2DShape(:,flippedIndexes(k));
    flippedColumn = flipud(bttb_col);
    BCCB_Column = [BCCB_Column; bttb_col; flippedColumn(2:N) ];
end
BCCB_ColumnsIn2DShape = reshape(BCCB_Column, [2*N, 2*N]);

IndexesBCCB = (1:2*N)';
Circulant_Indexes = [];
for i = 1:2*N
    Circulant_Indexes = [Circulant_Indexes, IndexesBCCB];
    IndexesBCCB = circshift(IndexesBCCB, 1);
end

BCCB_Matrix = [];
for l = 1:2*N
    ColumnIndexes = Circulant_Indexes(:,l);
    CirculantBlocks = [];
    for k = 1:2*N
        CurrentBlockIndex = ColumnIndexes(k);                                           % This is the block Index
        CirculantBlock = toeplitz(BCCB_ColumnsIn2DShape(:,CurrentBlockIndex));           % This is a circulant block now
        CirculantBlocks = [ CirculantBlocks ; CirculantBlock];                           % Append Column-wise C = [C1 ; C2; C3; C2]
    end
    BCCB_Matrix  = [BCCB_Matrix, CirculantBlocks];                                       % Append Row-wise
end
BCCB_Matrix

newY = BCCB_Matrix * sym_padded_x


%% MAking a reflection BTTB Matrix

REF_BTTB_Matrix = [BCCB_Matrix(3,:); BCCB_Matrix(4,:); BCCB_Matrix(1,:)]
