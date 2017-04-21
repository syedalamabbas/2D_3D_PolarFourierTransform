
%% This script is designed to test the adjoint operation of Spherical Polar Fourier Transform
% The adjoint of a transformation A is a unique transformation A^H so that
% <Ax,y> = <x, A^H y> for every x and y in the inner product space

%%
clear CheckPairsXX;
clear CheckPairsXZ;

clear
clc
addpath('..\NVIDIA_1DFrFT')
addpath('..\NVIDIA_2DPolarDFT')


%% Adjoint operation
N = 18;                                                % Always even
K = N+4;                                               % Always even
M = N+2;                                               % Always even

%% Impulse input 
% X = zeros(N+1, N+1, N+1)+ 1i*zeros(N+1, N+1, N+1);     % Complex matrix for computation of Spherical Polar Grid
% Y = zeros(K, M, N+1)  + 1i*zeros(K, M, N+1);           % Complex matrix of actual Spherical Polar Grid
% % Level 1 Check- This is the impulse response of the linear system
% X(1,1,1)= 1;
% Y(1,1,1) = 1;
% % Level 2 Check
%  X(1,1,:)= ones(1,N+1);
%  Y(1,1,:) = ones(1,N+1);
% %  Level 3 Check
% X(1,:,:)= randn(N+1,N+1);
% Y(1,:,:) = randn( M,N+1);
 
% X(:,13,:)= rand(N+1,N+1);
% Y(:,17,:) = rand(K,N+1);
 
% X(:,:,1)= ones(N+1,N+1);
% Y(:,:,1) = ones(K,M);

% % Level 4 with Real inputs only
% X(:,:,:)= randn(N+1,N+1,N+1);
% Y(:,:,:) = randn(K, M,N+1);


 

%% Final Random Complex Inputs
X = randn(N+1, N+1, N+1); %+ 1i*randn(N+1, N+1, N+1);                % Real Input matrix for computation of Spherical Polar Grid
Y = randn(K, M, N+1)+ 1i*randn(K, M, N+1);                        % Complex matrix of actual Spherical Polar Grid

% Direct computations  TEST
DAx = DirectBruteForce3D_DFT( X, K, M);
DAHy = DirectBruteForce3D_AdjointDFT( Y ); 

%% Checking Business for direct 
prod1  = squeeze(dot(DAx,Y));
prod2  = squeeze(dot(X,DAHy)); 
fprintf('\n The result of  (<Ax,y> - <x, A^H y>) comparison (direct computations) is :')
disp(abs(sum(sum(prod1)) -(sum(sum(prod2)))))
 
%% Vectorized version  TEST
Ax = VectorizedCompute3DSphericalDFT( X, K, M, [] );      % Output M x N+1 x N+1     Spherical Polar grid
AHy = Adjoint3DSphericalDFT(Y, []);                      % Output N+1 x N+1 x N+1   Adjoint image

disp('Max Difference between direct computations and my fast solution, forward and adjoint operators ');
disp(max(max(max(Ax-DAx))));
disp(max(max(max(AHy-DAHy))));
   
%% Checking Business for my vectorized
prod1  = squeeze(dot(Ax,Y)); 
prod2  = squeeze(dot(X,AHy));
fprintf('\n The result of  (<Ax,y> - <x, A^H y>) comparison (my solution) is :')
disp(abs(sum(sum(prod1)) -(sum(sum(prod2)))))
  