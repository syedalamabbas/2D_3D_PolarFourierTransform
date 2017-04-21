
%% This script is designed to test the adjoint operation of Cylindrical Transform
% The adjoint of a transformation A is a unique transformation A^H so that
% <Ax,y> = <x, A^H y> for every x and y in the inner product space

clear
clc
addpath('..\NVIDIA_2DPolarDFT')


%% Adjoint operation
N = 16;                                                % Always even
M = N+2;                                               % Always even
X = randn(N+1, N+1, N+1)+ 1i*randn(N+1, N+1, N+1);     % Complex matrix of Polar Grid
Y = randn(M, N+1, N+1)  + 1i*randn(M, N+1, N+1);       % Complex matrix
Ax = VectorizedCompute3DCylindricalDFT( X,  M );       % Output M x N+1 x N+1     Cylindrical grid
AHy = Adjoint3DCylindricalDFT(Y);                      % Output N+1 x N+1 x N+1   Adjoint image

%% Checking Business
prod1  = squeeze(dot(Ax,Y)) ;
prod2  = squeeze(dot(X,AHy));
fprintf('\n The result of  (<Ax,y> - <x, A^H y>) comparison is :')
disp(abs(sum(sum(prod1)) -(sum(sum(prod2)))))

   