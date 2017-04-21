clear all
clc


N = 88;   % always even
M = N+2;  % Always even again
X =  randn (N+1 );                            % Arbitrary Image Cartesian grid, this is a real image 
Y =  randn(M, N+1) + 1i* randn(M,  N+1);       % arbitrary Solution Polar grid , this must be complex values 
Ax = Compute2DPolarDFT (X, M);
AHy = Adjoint2DPolarDFT (Y);

prod1 = Y'.*Ax.';
prod2 = X'.*AHy;
disp(abs( sum(sum(prod1)) -conj(sum(sum(prod2)))));



X =  zeros (N+1 );   
X(1) = 1;         % Checking with impulse
[AxCorner,~,C] = Compute2DPolarCornersDFT (X, M);
YCorner =  randn(1, C) + 1i* randn(1, C);       % arbitrary values Polar grid corners, this must be complex values 
AHyCorner = Adjoint2DPolarCornersDFT (YCorner, N,M);

prod1Corner = YCorner'.*AxCorner.';  
prod2Corner = X'.*AHyCorner;
disp((abs( sum(sum(prod1Corner)) -conj(sum(sum(prod2Corner))))));




% X =  zeros (N+1 );   
% X(1) = 1;         % Checking with impulse
% [Ax,~,C] = Compute2DPolarCornersDFT (X, M);
% figure, imagesc(X)
% 
% % returnImage = Adjoint2DPolarCornersDFT (Ax, N,M);
% % figure, imagesc(real(returnImage),[0 255])
% %  
% Y =  randn(1, C) + 1i* randn(1, C);       % arbitrary values Polar grid corners, this must be complex values 
% AHy = Adjoint2DPolarCornersDFT (Y, N,M);
% 
% prod1 = Y'.*Ax.';  
% prod2 = X'.*AHy;
% % disp(eval(abs( sum(sum(prod1)) -conj(sum(sum(prod2)))))); % Use this for symbolic
% disp((abs( sum(sum(prod1)) -conj(sum(sum(prod2))))));