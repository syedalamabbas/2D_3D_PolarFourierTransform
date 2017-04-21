clear
clc
close all
addpath('..\NVIDIA_1DFrFT')
addpath('..\NVIDIA_2DPolarDFT')
addpath(genpath('..\..\..\Rotation Estimation'))
% addpath(genpath('C:\Users\Student\Documents\MATLAB\Dr.MahoorCode\Rotation Estimation'))

%% Load the original volume data
load('cat1_obj.mat');

% Before voxelization
PlotSurface1(vertices,faces);
title('Surface data ')
 
% After Voxelization
[bim] = verticestovolumefunc(vertices,faces);

origin=[0 0 0];
vxsize =[1 1 1];
[vertices, faces] =  gen_surf_data(bim,origin,vxsize);
PlotSurface1(vertices,faces);
title('Original Volume data')
% surfaceVal = 50;
surfaceVal = 1;
volumeOriginal = surfaceVal *double(bim);  % Fully Logical Array


%% Setup the inverse 
N = 54;  % Even
M = ( N+14);  % Even
K = (N+10);   % EEven 

% Image_3D = rand(N+1,N+1,N+1);
Image_3D = volumeOriginal;

SphericalPolarGrid  = VectorizedCompute3DSphericalDFT( Image_3D,K, M, []  ); 
Image_Final_double = Specialized3DInverse( SphericalPolarGrid );
 
Image_Final = (Image_Final_double > .56);       % Some adjustment 

[vertices, faces] =  gen_surf_data(Image_Final,origin,vxsize);
PlotSurface1(vertices,faces);
title('Reconstructed')
 
%% Error Plot 
error = (logical(Image_3D) - Image_Final);
disp('The max error is-')
disp(max(max(max(error))))
disp('The sum error is-')
disp(sum(sum(sum(error))))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the output
% The max error is-
%      1
% 
% The sum error is-
%      0
