

clear
clc
close all

%% Test Image data

TestImage = [ 0.1117 ,   0.5606 ,   0.6126 ,   0.8452  ,  0.8908 ,   0.4843 ,   0.4896   ,   0.5386 ,   0.2703,  2,2,;
    0.1363  ,  0.9296  ,  0.9900 ,   0.7386  ,  0.9823  ,  0.8449  ,  0.1925   , 0.6952  , 0.2085,  2,2,;
    0.6787  ,  0.6967  ,  0.5277  ,  0.5860 ,   0.7690 ,   0.2094  ,  0.1231   , 0.4991  ,  0.5650,  2,2,;
    0.4952  ,  0.5828  ,  0.4795  ,  0.2467  ,  0.5814  ,  0.5523  ,  0.2055 ,   0.5358  ,  0.6403, 2,2,;
    0.1897  , 0.8154  ,  0.8013  ,  0.6664  ,  0.9283  ,  0.6299  ,  0.1465  ,  0.4452   ,  0.4170, 2,2,;
    0.4950  ,  0.8790  ,  0.2278  ,  0.0835  ,  0.5801  ,  0.0320  ,  0.1891  ,  0.1239  ,  0.2060, 2,2,;
    0.1476  ,  0.9889  ,  0.4981  ,  0.6260  ,  0.0170  ,  0.6147  ,  0.0427  ,  0.4904  ,  0.9479, 2,2,;
    0.0550  ,  0.0005  ,  0.9009  ,  0.6609  ,  0.1209  ,  0.3624  ,  0.6352  ,  0.8530  ,  0.0821, 2,2,;
    0.8507  ,  0.8654  ,  0.5747  ,  0.7298  ,  0.8627  ,  0.0495  ,  0.2819  ,  0.8739  ,  0.1057 , 2,2;
    0.5386, 0.2703, 0.1420, 0.8604, 0.0309, 0.5590, 0.1182, 0.8182, 0.9052, 2, 2;
    0.6952, 0.2085, 0.1665, 0.9344, 0.9391, 0.8541, 0.9884, 0.1002, 0.6754, 2, 2
    ];



%% Defining Image

Size = 251;            % Odd
% I = double(phantom('Modified Shepp-Logan',Size));
% I = double(imresize(imread('rice.png'), [Size,Size]));
%  I = double(imresize(imread('cameraman.tif'), [Size,Size]));
%  I = double(imresize(imread('mandi.tif'), [Size,Size]));
 I = double(imresize(imread('lena.tif'), [Size,Size]));
%  I = randn (Size,Size); 


%% MATLAB Filterred Back Projection
P =  I; 
THETA = 0:.7:179;
R = radon(P,THETA);
I1 = iradon(R,THETA);
I2 = iradon(R,THETA,'linear','none');
figure, imshow(P, []), title('Original')
figure, imshow(I1,[]), title('Filtered backprojection')
figure, imshow(I2,[]), title('Unfiltered backprojection')
figure, imagesc((P-imresize(I1, size(P)))) , title('Error with FBP')
%% My algorithm
% padVal = 91;
% I = padarray(I, [padVal padVal], 'both');

[Size,~] = size(I);

figure , imshow(I, [0 255]);
% figure , imagesc(I);
% colorbar 
title ('Original Image')
% I = randn (Size,Size); 
% I = TestImage;
N = Size -1;      % N is even

  

%% Testing forward and reverse FrFT on the butterfly grid
% gridSpacing =  -N/2:N/2;
% newSpacing = abs(gridSpacing)* sind(180/M);
% newSpacing (N/2 +1) = 1;  
% alpha = cosd(180/M);
 

%% 1D FrFT test forward and reverse
% vector = rand(21,1 );
% FrFTForwardCol = VectorizedFrFT_Centered(vector,alpha);
% IFrFTForwardCol = VectorizedFrFT_Centered(FrFTForwardCol,-alpha)/ length(vector);
%
% figure, plot(1:length(vector),vector, 1:length(vector), (real(IFrFTForwardCol)))
% legend('Original', 'Reconstructed')
% grid on
%
% figure, plot (1:length(vector),vector - (real(IFrFTForwardCol)))
% legend('Error 1D FrFT')
% grid on
%
% FFTForwardCol = VectorizedFrFT_Centered(vector,1);
% IFFTForwardCol = VectorizedFrFT_Centered(FrFTForwardCol,-1)/ length(vector);
%
% figure, plot (1:length(vector),vector -(real(IFrFTForwardCol)))
% legend('Error 1D FFT using FrFT code')
% grid on
%
%  figure , plot (1:length(vector),vector - ifft(fft(vector)))
% legend('Error 1D FFT')
% grid on
%
% %% 2D FrFT test
% simpleFFT = fftshift(fft2(TestImage));
% forwardFrFTButterfly = ReconstructFromButterFlyData( TestImage,  -newSpacing, -alpha);
%
% simpleReverseFFT = (ifft2(ifftshift(simpleFFT)));
% InverseFrFTButterfly = ReconstructFromButterFlyData( forwardFrFTButterfly,  newSpacing, alpha)/(N+1)^2;
%
% figure, imagesc(real(InverseFrFTButterfly));
% xlabel('x')
% ylabel('y');
% title('Reconstructed Image')




%% Testing the inverse polar 2D DFT


noOfAngles = 2*(Size +1);
Polar_Grid = VectorizedCompute2DPolarDFT( I,  noOfAngles ) ;
[Polar_corners,~,C]  = Compute2DPolarCornersDFT2( I,  noOfAngles );
[ ReconstImage ] =  Inverse2DPolarDFT( Polar_Grid, Polar_corners); 

% [ ReconstImage ] =  Inverse2DPolarWithCorners( Polar_Grid, Polar_corners);  % output is real
 
% [ ReconstImage ] = DirectInverse2DPolarDFT( Polar_Grid );  % output is real
    
 figure,     
 hold on
subplot(2,2,1) % first subplot     
axis ([-Size Size -Size Size]) 
imshow(I, [0 255]);
% imshow(I)
xlabel('x')
ylabel('y');  
title('Original Image')  
subplot(2,2,2) % second subplot
axis equal
imshow(real(ReconstImage), []);  
xlabel('x') 
ylabel('y');
title('Reconstructed Image fast exact solution')
disp('The difference norm between the exact and reconstructed is ')
disp(norm(real(ReconstImage)-I)); 
  
figure, imagesc(real(ReconstImage)-I) ; colorbar
title('Reconstruction Error as fixed pattern noise')


%% Testing generic matlab functionality
% myN =  Size^2;
% A = complex(rand(myN,myN),rand(myN,myN)) ; %gallery('wilk',21);
% b = complex(rand(myN,1),rand(myN,1)) ; sum(A,2);
% gmres_x = gmres(A,b, myN,10^-6,myN);
% cg_x = cgs(A'*A,A'*b, 10^-6,myN);
    
    
%% Checking the pseudo polar behaviour  
      
X= imresize(I, [N,N]); %randn(N);  +sqrt(-1)*randn(N);   
Y=PPFFT(X,1,1); 
X1=IPPFFT(Y);     
   
% figure, 
subplot(2,2,3) % first subplot
imshow(X, [0 255]); 
xlabel('x')
ylabel('y');
title('Original Image')
subplot(2,2,4) % second subplot
axis equal
imshow(real(X1), [0 255]);
xlabel('x')
ylabel('y');
title('Reconstructed Image using PP-grid')
hold on
disp('The difference norm between the exact and reconstructed is ')
disp(norm(X1-X))


figure, 
subplot(1,2,1)
imagesc(real(ReconstImage)-I)
xlabel('x')
ylabel('y');
title('Reconstruction Error my solution')
subplot(1,2,2)
axis equal
imagesc(real(X1)-X)
xlabel('x')
ylabel('y');
title('Reconstruction Error PP grid')




% %% Comparison
%
% theta = 23.46;   % angular location in degrees
% r = 3;           % picking the 3rd circle
% desiredPoint = [r*cosd(theta) r*sind(theta) ];
%
% DFT_Point = DirectSinglePoint_2DDFT( I, desiredPoint )   % Direct computation
%
% Fourier_Coefficients_Circle_half = Polar_Grid (:,N/2+1-r );
% Fourier_Coefficients_Circle_otherhalf = Polar_Grid (:,N/2+1+r );
% Fourier_Coefficients_Circle = [(Fourier_Coefficients_Circle_half)'  (Fourier_Coefficients_Circle_otherhalf)' ];
%
% Interp_Complex_Point = FourierBasedTrigInterp( desiredPoint , fliplr(Fourier_Coefficients_Circle)  , noOfAngles )