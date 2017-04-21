
 clc
 clear 
 close all
N = 3;
a = sym('a_', [1 N^2]);



ImageFinal = zeros(N,N);
ImageFinal(1,1) = 1;            % Making it an impulse or a point spread function, 
BTTB_ColumnsIn2DShape = Adjoint2DPolarDFT( Compute2DPolarDFT( ImageFinal,  N+1 ));     % It is (N+1) x (N+1) size, it is a impluse response of the linear system

BTTB_ColumnsIn2DShape = real(BTTB_ColumnsIn2DShape);
Indexes = 1:(N);
Toeplitz_Indexes = toeplitz(Indexes);   % Block Level

BTTB_Matrix = []; 
for l = 1:(N)
    ColumnIndexes = Toeplitz_Indexes(:,l);
    ToeplitzBlocks = [];
    for k = 1:(N)
        CurrentBlockIndex = ColumnIndexes(k);                      % This is the block Index
        ToeplitzBlock = toeplitz(BTTB_ColumnsIn2DShape(:,CurrentBlockIndex));
        ToeplitzBlocks = [ ToeplitzBlocks ; ToeplitzBlock];    % Append Column-wise T = [T1 ; T2; T3]
    end
    BTTB_Matrix  = [BTTB_Matrix, ToeplitzBlocks];              % Append Row-wise
end

figure, imagesc(BTTB_Matrix)
title('BTTB Matrix')
%% prodReal is a simplified product computed using script Symbolic Math for Butterfly grid , so we must run it first
NumericProdReal = BTTB_Matrix;
decimalPlaces = 1000000; 
% BTTB  
Row = NumericProdReal (1,:);
Row = round((Row*decimalPlaces))/decimalPlaces;
for m=1: N^2 
    for n = 1: N^2
%         m
%         n
           disp('Just finished iteration, (m,n) ')
           disp(m)
           disp(n)
         currentElement = NumericProdReal(m,n);
         currentElement = round((currentElement*decimalPlaces))/decimalPlaces;
         indexes = find (Row  == currentElement );
         BTTB(m,n) = a(indexes(1)); 
    end
end 

disp('This is the Block Toeplitz Toeplitz Block Matrix !')
 BTTB 
 
% Toeplitz individual blocks

T = sym('T', [1 N]);
counter = 0;
disp('These are the individual Toeplitz Blocks Matrix !')
for k = 1:N
    offset = counter*N ;
    disp(T(k) )
    disp(BTTB(1:N, (1+ offset):(N+offset)) )
    counter = counter +1;
end


%% For 25 x 25 BTTB
if (N == 5)
syms T1 T2 T3 T4 T5
Composite = [ T1, T2, T3, T4, T5 ; T2 , T1, T2, T3, T4; T3, T2, T1, T2, T3; T4, T3, T2, T1, T2  ; T5, T4, T3, T2, T1 ]
end

if (N == 3)
    syms T1 T2 T3
    Composite = [ T1, T2, T3; T2 , T1, T2; T3, T2, T1]
end


% syms a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20 a21 a22 a23 a24 a25
% 
% T1 = [ a1, a2, a3, a4, a5;
%  a2, a1, a2, a3, a4;
%  a3, a2, a1, a2, a3;
%  a4, a3, a2, a1, a2;
%  a5, a4, a3, a2, a1];
% T2 = [  a2, a7, a8, a9, a10;
%   a7, a2, a7, a8,  a9;
%   a8, a7, a2, a7,  a8; 
%   a9, a8, a7, a2,  a7;
%  a10, a9, a8, a7,  a2];
% T3 = [  a3,  a8, a13, a14, a15;
%   a8,  a3,  a8, a13, a14;
%  a13,  a8,  a3,  a8, a13;
%  a14, a13,  a8,  a3,  a8;
%  a15, a14, a13,  a8,  a3];
% T4 = [  a4,  a9, a14, a19, a20;
%   a9,  a4,  a9, a14, a19;
%  a14,  a9,  a4,  a9, a14;
%  a19, a14,  a9,  a4,  a9;
%  a20, a19, a14,  a9,  a4];
% T5 =  [  a5, a10, a15, a20, a25;
%  a10,  a5, a10, a15, a20;
%  a15, a10,  a5, a10, a15;
%  a20, a15, a10,  a5, a10;
%  a25, a20, a15, a10,  a5];

% disp('Difference between composite and BTTB') 
% disp(eval(Composite)-BTTB)        % Should be all Zeros if the arrangement is correct !


%% Verify if these matrices are indeed Toeplitz
% T1 - toeplitz (T1(1,:))
% T2 - toeplitz (T2(1,:))
% T3 - toeplitz (T3(1,:))
% T4 - toeplitz (T4(1,:))
% T5 - toeplitz (T5(1,:))  % Should be all Zeros

%% Generating the Block Circulant Matrix
if (N == 3)
    ColCirc  = [Composite(:,1); T2];
end
if (N == 5)
    ColCirc  = [Composite(:,1); T4 ; T3 ; T2];
end

Circulant = [];
for i = 1:length(ColCirc)
    Circulant = [Circulant, ColCirc];
    ColCirc = circshift(ColCirc, 1);
end
disp('Block Circulant=')
disp(Circulant)

%% Generating Block Circulant Circular Blocks
if (N == 3)
 syms ref_T1 ref_T2 ref_T3

syms a_1 a_2 a_3

T1 =[ a_1, a_2, a_3;
 a_2, a_1, a_2;
 a_3, a_2, a_1];

ref_T1 = [   0, a_3, a_2;
 a_3,   0, a_3;
 a_2, a_3,   0];

CirculantBlock =[  T1, ref_T1; ref_T1, T1]
 

BCCB = [ T1, ref_T1, T2, ref_T2, T3, ref_T3, T2, ref_T2;
  ref_T1, T1, ref_T2, T2, ref_T3, T3,  ref_T2,  T2 ;
  T2, ref_T2, T1, ref_T1, T2, ref_T2, T3, ref_T3;
  ref_T2, T2,  ref_T1, T1,  ref_T2, T2,  ref_T3, T3;
  T3,  ref_T3, T2, ref_T2, T1, ref_T1, T2, ref_T2;
  ref_T3, T3,   ref_T2, T2,  ref_T1, T1,  ref_T2, T2 ;
  T2, ref_T2, T3,  ref_T3, T2, ref_T2, T1, ref_T1;
  ref_T2, T2,  ref_T3, T3,   ref_T2, T2,  ref_T1, T1  ]
end

