close all
clc

%% Initializations
N = 18;    % Always even
M= 2*(N+1);
I = randn(N+1);   % Input Image

figure , imagesc(I);
title ('Original Image')
colorbar

% %% Padding the image
% NewN_2 = floor(N/2/(cosd(45)));
% NewSize = 2*NewN_2+1;
% 
% I_New =  padarray (I, [NewN_2-N/2, NewN_2-N/2], 'both' );
% figure , imagesc(I_New);
% title ('Padded Image')
% 
% % Rehashing the solution now
% I = I_New;
% N = 2*NewN_2;
% M = 2*(N+1);

%% Computing the Forward Transform
Polarx2_withCartesianCornerTransFormMatrix =  CreateMatrixFromGrid( 'EP_CartC', N ,M,0);    % This is Polar x2 radially with cartesian corners
I_1D = reshape(I,[(N+1)^2,1]);
DFTValues = Polarx2_withCartesianCornerTransFormMatrix * I_1D;


%% Nulling the corner values !!
% [sizeX,sizeY ] = size(Polarx2_withCartesianCornerTransFormMatrix);
% lastIndexPolar = M*(2*N+1);
% array = lastIndexPolar+1:sizeX;
% DFTValues(array)= zeros([length(array),1]) ;

%% Computing inverse.

A_p = Polarx2_withCartesianCornerTransFormMatrix.'*Polarx2_withCartesianCornerTransFormMatrix; % A^H A
figure, imagesc(real(A_p))
title('Two level Toeplitz matrix')
colorbar

SVDs = svd(real(A_p));
figure, plot(1:length(SVDs),flipud(SVDs))
title(strcat('The singular values decay, the condition number of A_p = ',num2str(cond(A_p),4)))

b_p = Polarx2_withCartesianCornerTransFormMatrix.'* DFTValues;                                 % A^H b   for Ax = b system

ImageFinal_1D = cgs(A_p , b_p, [], 200);
ImageFinal = reshape(ImageFinal_1D,[(N+1),(N+1)]);
figure, imagesc(I - real(ImageFinal))
title('Error in reconstruction')
colorbar
figure,imagesc(real(ImageFinal))
title('Reconstructed image')
colorbar


%% Testing the two structure like A = A_1 : A_2 , A^H * A = A^H_1 * A_1 + A^H_2 * A_2
N = 8;
a= sym('a',[1,N]);
b= sym('b',[1,N]); 
c= sym('c',[1,N]); 
d= sym('d',[1,N]);
e= sym('e',[1,N]); 
f= sym('f',[1,N]);
g= sym('g',[1,N]);
h= sym('h',[1,N]);
i = sym('i',[1,N]);
j = sym('j',[1,N]);

A_1 = [a;b;c;d;e]
A_2 = [f;g;h;i;j]
A = [a;b;c;d;e;f;g;h;i;j]
A'*A - (A_1'*A_1 + A_2'*A_2)   % The error to be zero
 
