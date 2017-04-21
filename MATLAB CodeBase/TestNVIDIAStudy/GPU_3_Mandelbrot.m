clear;
clc;

close all;
% A 1000x1000 grid of real parts (X) and imaginary parts (Y) is created 
% between these limits and the Mandelbrot algorithm is iterated at each grid location. 
% For this particular location 500 iterations will be enough to fully render the image.
maxIterations = 500;
gridSize = 1200;
xlim = [-0.748766713922161, -0.748766707771757];
ylim = [ 0.123640844894862,  0.123640851045266];

%% Direct and Simple Mandelbrot Set in MATLAB using CPU only

% Setup
t = tic();
x = linspace( xlim(1), xlim(2), gridSize );
y = linspace( ylim(1), ylim(2), gridSize );
[xGrid,yGrid] = meshgrid( x, y );
z0 = xGrid + 1i*yGrid;
count = ones( size(z0) );

% Calculate
z = z0;
for n = 0:maxIterations
    z = z.*z + z0;
    inside = abs( z )<=2;
    count = count + inside;
end
count = log( count );
  
% Show 
cpuTime = toc( t );  
fig = gcf; 
fig.Position = [200 200 600 600]; 
imagesc( x, y, count );
axis image
colormap( [jet();flipud( jet() );0 0 0] );
title( sprintf( '%1.2fsecs (without GPU)', cpuTime ) );

%% Using GPUArray

% Setup 
t = tic();
x = gpuArray.linspace( xlim(1), xlim(2), gridSize );
y = gpuArray.linspace( ylim(1), ylim(2), gridSize );
[xGrid,yGrid] = meshgrid( x, y );
z0 = complex( xGrid, yGrid );
count = ones( size(z0), 'gpuArray' ); 
 
% Calculate
z = z0;
for n = 0:maxIterations
    z = z.*z + z0;
    inside = abs( z )<=2;
    count = count + inside;
end
count = log( count );

% Show
count = gather( count ); % Fetch the data back from the GPU
naiveGPUTime = toc( t );
figure,
fig = gcf;
fig.Position = [200 200 600 600];
imagesc( x, y, count )
axis image
colormap( [jet();flipud( jet() );0 0 0] );
title( sprintf( '%1.3fsecs (naive GPU) = %1.1fx faster', ...
    naiveGPUTime, cpuTime/naiveGPUTime ) )

%% Using Arrayfun

% Setup
t = tic();
x = gpuArray.linspace( xlim(1), xlim(2), gridSize );
y = gpuArray.linspace( ylim(1), ylim(2), gridSize );
[xGrid,yGrid] = meshgrid( x, y );

% Calculate
count = arrayfun( @pctdemo_processMandelbrotElement, ...
                  xGrid, yGrid, maxIterations );

% Show
count = gather( count ); % Fetch the data back from the GPU
gpuArrayfunTime = toc( t );
figure,
fig = gcf;
fig.Position = [200 200 600 600];
imagesc( x, y, count ) 
axis image
colormap( [jet();flipud( jet() );0 0 0] );
title( sprintf( '%1.3fsecs (GPU arrayfun) = %1.1fx faster', ...
    gpuArrayfunTime, cpuTime/gpuArrayfunTime ) );


%% Using ArrayFire GPU Accelerated Library
% addpath('D:\GitHub\arrayfire-windows-scripts\SimpleCUDAProj\MandelbrotCUDATest1\x64\Debug')
% addpath('D:\GitHub\MexCUDA_ArrayFire\x64\Release')
addpath('D:\GitHub\MexCUDA_ArrayFire\x64\Debug')
t = tic();
x = gpuArray.linspace( xlim(1), xlim(2), gridSize );
y = gpuArray.linspace( ylim(1), ylim(2), gridSize );
[xGrid,] = meshgrid( x, y );
count = ones( size(xGrid), 'gpuArray' );
Outcount = MEXMandelbrotCUDAArrayFire(x,y, count,maxIterations, gridSize);
% count = gather( count ); % Fetch the data back from the GPU
gpuArrayFireTime = toc( t );  
figure, 
fig = gcf;   
fig.Position = [200 200 600 600]; 
imagesc( x, y, Outcount' ) 
axis image 
colormap( [jet();flipud( jet() );0 0 0] );
title( sprintf( '%1.3fsecs (With ArrayFire GPU Accelerated Library) = %1.1fx faster', ...
    gpuArrayFireTime, cpuTime/gpuArrayFireTime ) );

%%  Using CUDA kernel 

% Load the kernel 
cudaFilename = 'mandelbrotViewerProcessElement.cu';
ptxFilename =  ['mandelbrotViewerProcessElement.',parallel.gpu.ptxext];
kernel = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename );

% Setup 
t = tic();
x = gpuArray.linspace( xlim(1), xlim(2), gridSize );
y = gpuArray.linspace( ylim(1), ylim(2), gridSize );
[xGrid,yGrid] = meshgrid( x, y );

% Make sure we have sufficient blocks to cover all of the locations
numElements = numel( xGrid );
kernel.ThreadBlockSize = [kernel.MaxThreadsPerBlock,1,1];
kernel.GridSize = [ceil(numElements/kernel.MaxThreadsPerBlock),1];

% Call the kernel
count = zeros( size(xGrid), 'gpuArray' );
count = feval( kernel, count, xGrid, yGrid, 4.0,maxIterations, numElements );

% Show
count = gather( count ); % Fetch the data back from the GPU
gpuCUDAKernelTime = toc( t );
figure,
fig = gcf;
fig.Position = [200 200 600 600];
imagesc( x, y, count )
axis image
colormap( [jet();flipud( jet() );0 0 0] );
title( sprintf( '%1.3fsecs (GPU CUDAKernel) = %1.1fx faster', ...
    gpuCUDAKernelTime, cpuTime/gpuCUDAKernelTime ) );