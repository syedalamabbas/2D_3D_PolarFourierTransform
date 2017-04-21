

clear
clc
close all

addpath('..\NVIDIA_1DFrFT')

%% Preparing the image
N = 250;  % Always even
% M = 2*(N+1); % Number of angles always even too
I = double(imresize(imread('lena.tif'), [N+1,N+1]));
% I = double(imresize(imread('cameraman.tif'), [N+1,N+1]));
% I = randn(N+1);
figure, imshow(I, [0 255])
title('Original Image')

NewN_2 = floor(N/2/(cosd(45)));
NewSize = 2*NewN_2+1;

I_New =  padarray (I, [NewN_2-N/2, NewN_2-N/2], 'both' );
figure , imshow(I_New, [0 255]);
title ('Padded Image')


figure,  imshow (I_New(NewN_2-N/2:NewN_2+N/2,NewN_2-N/2:NewN_2+N/2 ), [0 255] )
title('UnPadded ! ')



% Restructuring 
I = I_New;
N = 2*NewN_2;
M = 2*(N+1);
M = 30;     % Sparse view ??? Need better regularization

%% Forward  construction  same as full grid but retaining only Corners
I_DFT = fftshift(fft2(I));
figure, imagesc(log(1+abs(I_DFT)))
title('Image, magnitude of Fourier Transform')
colormap hot

%% Polar Nulling
newGridSpacing = -N/2:N/2;
[CartesianGridX, CartesianGridY]= meshgrid(newGridSpacing,newGridSpacing);
Weight=(CartesianGridX.^2 + CartesianGridY.^2 >= (N/2)^2);
Pos=find(Weight(:));
zeroArray = zeros(length(Pos),1);
I_DFT(Pos) = zeroArray;
figure, imagesc(log(1+abs(I_DFT)))
% title('Image, magnitude of Fourier Transform after the central polar nulling ')
colormap hot
% axis equal
axis tight
 I_DFT = zeros(N+1);  

%% Proper inverse
PolarGrid  = Compute2DPolarDFT( I ,  M ); 
Image_Final = Specialized2DInverse( PolarGrid,I_DFT, N, M );

figure, imshow(Image_Final, [0 255]) 
title('Reconstructed Image (Exact)')
 
%% Adjoint or Inverse reconstruction
% I_Reconstruction = ifft2(ifftshift(I_DFT));
% figure, imagesc(I_Reconstruction) 
% title('Reconstructed Image') 

error = Image_Final - I;
figure, imagesc(error )
title('Error (Exact)')
% colormap cool



%% Testing with MATLAB radon and iradon function for forward and inverse solution

P =  I; 
THETA = 0:.35:179;
R = radon(P,THETA); 
I1 = iradon(R,THETA);
I2 = iradon(R,THETA,'linear','none');
figure, imshow(P, []), title('Original')
figure, imshow(I1,[]), title('Reconstructed Image (FBP)')
figure, imshow(I2,[]), title('Unfiltered backprojection')
figure, imagesc((P-imresize(I1, size(P)))) , title('Error (FBP)')

 




%% Adding now polar Compute the toeplitz BTTB
% ImageFinal = zeros(N+1);
% I_DFT = fftshift(fft2(ImageFinal));
% I_DFT(Pos) = zeroArray;
% ComputePolar_CartesianCorner = @(X) ;
% AH_A =  Adjoint2DPolarDFT((Compute2DPolarDFT( ImageFinal,  M ) )) + ifft2(ifftshift(I_DFT));


