clear;
clc;
close all

addpath('..\NVIDIA_1DFrFT')
addpath('..\NVIDIA_2DPolarDFT')
addpath('..\MEXArrayFireCUDA')

%% Large data with 45 degrees angle included
Size = 19;   % Odd
I = double(rand(Size,Size, Size));
% I = get3D_9_by_9StoredData();
[sizeX,sizeY,sizeZ] =  size(I);
N = sizeX -1;      % N is even
% I(:,:,1)
noOfAnglesTheta = Size+1;       % Always even
noOfAnglesPhi = Size+1;         % Always even
 
%% Direct Computation brute force
disp('Testing Brute Force vs My Fast Algorithm ')
tic
LHS = DirectBruteForce3D_DFT( I,noOfAnglesTheta,noOfAnglesPhi);
disp(['Execution time on CPU (direct) = ',num2str(toc), ' (in seconds)']);

%% Fast Spherical transform 
M = noOfAnglesPhi;          % M is also even
K = noOfAnglesTheta;       % K is also even
tic  
RHS_CPU = VectorizedCompute3DSphericalDFT(I, K, M,[]);
disp(['Execution time on CPU (my special spherical fast) = ',num2str(toc), ' (in seconds)']);
disp( ['Maximum absolute error sum = ',num2str(max(max(max(abs(sum(LHS - RHS_CPU))))))])
 
%% GPU   
tic   
RHS_GPU = VectorizedCompute3DSphericalDFT(gpuArray(I), K, M,[]);
f = @()VectorizedCompute3DSphericalDFT(gpuArray(I), K, M,[]);

gpuDirectTime  = toc;
disp(['Execution time on GPU (my special spherical fast) = ',num2str(gpuDirectTime), ' (in seconds)']);
disp( ['Maximum absolute error sum = ',num2str(max(max(max(abs(sum(LHS - RHS_GPU))))))]);
gpuDirectTime = gputimeit(f);
disp(['Execution time using gpuTimeIT = ', num2str(gpuDirectTime), ' (in seconds)'])
 
%% Specialized MEX_CUDA_ArrayFire Function call

%  Using arrayfire mex, capable on 3 platforms, currently designed to use only the CUDA 
CURRENT_PLAT = PLATFORMS.CUDA;       % This is zero by default
 
switch(CURRENT_PLAT) 
    case PLATFORMS.CUDA % CUDA Path addition
%          addpath('D:\GitHub\RadonTransform\SphericalPolar3DTransformMEXArrayFire\x64\Release\');
%          addpath('D:\GitHub\RadonTransform\SphericalPolar3DTransformMEXArrayFire\x64\Release\');
         addpath('..\MEXArrayFireCUDA');
         gI = gpuArray(I);                               % Move array to GPU
    case PLATFORMS.CPU % CPU path addition
%         addpath('D:\GitHub\RadonTransform\TwoDPolarTransformMEXArrayFire\x64\CPU_Debug');
%         addpath('D:\GitHub\RadonTransform\TwoDPolarTransformMEXArrayFire\x64\CPU_Release');
%         gI = I;
    case PLATFORMS.OPENCL
end 
 
t = tic();
[realImage, ImagImage]= SphericalPolar3DTransformMEXArrayFire(gI, K, M,Size, uint32( CURRENT_PLAT));
mexOutput = complex(realImage,ImagImage);

g1time =  toc( t );
if(CURRENT_PLAT == PLATFORMS.CUDA)
    clear mex;  
end  
fprintf('Execution time on GPU (ArrayFire, fully vectorized! on CUDA) = %1.3f secs, it is %1.1fx faster ! \n', ...
    g1time, gpuDirectTime/g1time ); 
disp( ['Maximum absolute error sum = ',num2str(max(max(max(abs(sum(RHS_GPU - mexOutput))))))]);
disp(['Execution time using gpuTimeIT = ', num2str(g1time), ' (in seconds)'])


%% Direct Computation brute force with single precision
tic
LHS_Single = SinglePrecDirectBruteForce3D_DFT( I,noOfAnglesTheta,noOfAnglesPhi);
disp(['Execution time on CPU (direct) = ',num2str(toc), ' (in seconds)']);

%% Specialized MEX_CUDA_ArrayFire Function call with single precision

%  Using arrayfire mex, capable on 3 platforms, currently designed to use only the CUDA 
CURRENT_PLAT = PLATFORMS.CUDA;       % This is zero by default
 
switch(CURRENT_PLAT) 
    case PLATFORMS.CUDA % CUDA Path addition
%          addpath('D:\GitHub\RadonTransform\SphericalPolar3DTransformMEXArrayFire\x64\Release\');
%          addpath('D:\GitHub\RadonTransform\SphericalPolar3DTransformMEXArrayFire\x64\Release\');
         addpath('..\MEXArrayFireCUDA');
         gI = gpuArray(single(I));                               % Move array to GPU
    case PLATFORMS.CPU % CPU path addition
%         addpath('D:\GitHub\RadonTransform\TwoDPolarTransformMEXArrayFire\x64\CPU_Debug');
%         addpath('D:\GitHub\RadonTransform\TwoDPolarTransformMEXArrayFire\x64\CPU_Release');
%         gI = I;
    case PLATFORMS.OPENCL
end 
 
t = tic();
[realImage, ImagImage]= SphericalPolar3DTransformMEXArrayFireFloatPrec(gI, K, M,Size, uint32( CURRENT_PLAT));
mexOutput = single(gather(complex(realImage,ImagImage)));

g1time =  toc( t );
if(CURRENT_PLAT == PLATFORMS.CUDA)
    clear mex;  
end  
% fprintf('Execution time on GPU (ArrayFire, fully vectorized! on CUDA single precision) = %1.3f secs, it is %1.1fx faster ! \n', ...
%     g1time, gpuDirectTime/g1time ); 
disp( ['Maximum absolute error  = ',num2str(max(max(max(abs((LHS_Single - mexOutput))))))]);
disp( ['Maximum absolute error sum = ',num2str(max(max(max(abs(sum(LHS_Single - mexOutput))))))]);
disp(['Execution time using gpuTimeIT = ', num2str(g1time), ' (in seconds)'])
 
%% Capture CPU / GPU performance
Sizes = 15:2:168; 
TotalNoOfSizes = length(Sizes);
runningTimeCPU = zeros(1,TotalNoOfSizes );
runningTimeGPU = zeros(1,TotalNoOfSizes );

for k = 1:TotalNoOfSizes
    disp(['Iteration',num2str(k)]);
    
    size = Sizes(k);
    M = size + 29;     % M is also even
    K = size + 13;       % K is also even
    I = double(rand(size,size, size));
    
    tic
    RHS_CPU = VectorizedCompute3DSphericalDFT(I, K, M,[]);
    runningTimeCPU(k) = toc;
    
    f = @()VectorizedCompute3DSphericalDFT(gpuArray(I), K, M,[]);
    runningTimeGPU(k) = gputimeit(f);
end
csvwrite('runningTimeCPU.dat',runningTimeCPU)
csvwrite('runningTimeGPU.dat',runningTimeGPU)

type runningTimeCPU.dat
type runningTimeGPU.dat