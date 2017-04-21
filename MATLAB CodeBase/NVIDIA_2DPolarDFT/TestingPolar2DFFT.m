
clear;
clc;
%
% I = rand (9,9);
% [sizeX,sizeY] =  size(I);
% N = sizeX -1;      % N is even
% F_x_alpha_i = zeros(size(I));   % X- scaled square grid
% F_y_alpha_i = zeros(size(I));   % Y- scaled square grid
%
% % X axis Scaling
% for x = 1:N+1
%     row = I(x,:);
%     FrFtSeq = FrFT_Centered(row, 1.2);
% %     FrFtSeq = fftshift(fft(row));
%     F_x_alpha_i(x,:) = FrFtSeq;
% end
%
% % Y axis Scaling
% for y = 1:N+1
%     col = F_x_alpha_i(:,y);
%     FrFtSeq = FrFT_Centered(col, .4);
% %      FrFtSeq = fftshift(fft(col));
%     F_y_alpha_i(:,y) = FrFtSeq;
% end
%
% K = fftshift(fft2(I));
%
% disp('Testing my FrFTCentered code again, just the difference in the DC value with Direct Computation using FFT2')
% F_y_alpha_i (N/2+1, N/2+1) - K(N/2+1, N/2+1)
%
%
% disp('Testing Brute Force vs My Fast Algorithm ')
% LHS = DirectBruteForce2D_DFT( I, N+2);
% RHS = Compute2DPolarDFT( I,  N+2);
% max(abs(LHS-RHS))
%
% LHS
%
% RHS

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


%% Large data with 45 degrees angle included
%  Size = 1601;   % Odd

% Size = 69;
Size = 21;            % Odd
 I = double(rand (Size,Size));
%  I = TestImage;
[sizeX,sizeY] =  size(I);
N = sizeX -1;      % N is even

addpath('..\NVIDIA_1DFrFT')

%% Direct Computations
disp('Testing Brute Force vs My Fast Algorithm ')
LHS = DirectBruteForce2D_DFT( I, Size+1);
ctime = timeit(@()VectorizedCompute2DPolarDFT( I,  Size+1));       % Measure CPU time
disp(['Execution time on CPU (direct) = ',num2str(ctime), ' secs.']);
 

%% Fast Non-vectorized on CPU
RHS = Compute2DPolarDFT( double(I),  Size+1);
c1time = timeit(@()Compute2DPolarDFT( I,  Size+1));       % Measure CPU time
fprintf('Execution time on CPU (myFast, non-vector) = %1.3f secs, it is %1.1fx faster ! \n', ...
    c1time, ctime/c1time );
disp(['Maximum absolute errors element wise = ',num2str(max(max(abs(LHS-RHS))))]);
  
     
%% Fast Vectorized on CPU
 NewRHS = VectorizedCompute2DPolarDFT( I,  Size+1);
c2time = timeit(@()VectorizedCompute2DPolarDFT( I,  Size+1));       % Measure CPU time
fprintf('Execution time on CPU (myFast, half way vectorized) = %1.3f secs, it is %1.1fx faster ! \n', ...
    c2time, ctime/c2time );
disp(['Maximum absolute errors element wise  = ',num2str(max(max(abs(LHS-NewRHS))))]);
% disp(['Maximum relative errors element wise  = ',num2str(max(max(abs(LHS-NewRHS)))/ max(max(abs(LHS))))]);
%  


global BetaMaps; 
global E_ns;
global PremultiplicationFactors;
global PostmultiplicationFactors;
global Z ;

%% Vectorized and Precomputations on CPU


t = tic();
gI = I;
M = Size+1;                                    % Number of angles
L = (M-2)/4;                                   % Number of levels
hasFortyFiveDegrees = 0;
if(rem(M-2,4) ~= 0)
    hasFortyFiveDegrees = 1;
    L = ceil (L);                % use + 1 to compute for 45
end

n=[0:1:N, -N-1:1:-1]';
J = (0:1:N)';
K = (0:1:N)';
Z = zeros(N+1,N+1, 'like',I);
J = J - N/2;

alphas = cosd((1:1:L)*180/M);
betas = sind((1:1:L)*180/M);

LineSpacing = (-N/2:N/2);

BetaMaps = zeros ( N+1,N+1, 2*L,'like',I);
E_ns = zeros (L, 2*N+2,'like',I);
PremultiplicationFactors = zeros (L+1, N+1,'like',I);
PostmultiplicationFactors = zeros (L+1, N+1,'like',I);

lastIndex = 0;
for l=1:L
    E_ns (l,:) =  exp(-1i*pi*alphas(l)*n.^2/(N+1));                  % Sequence as defined in the paper
    PremultiplicationFactors(l,:) = exp(1i*pi*J *N*alphas(l)/(N+1));
    PostmultiplicationFactors(l,:) = exp(1i*pi*alphas(l)*N*K/(N+1));
    firstLocations =  (1 + lastIndex): lastIndex + (N+1);
    conjLocations = (1 + lastIndex+ (N+1)): lastIndex + 2*(N+1);
    lastIndex = lastIndex + 2*(N+1);
    BetaMaps( :,firstLocations ) =  exp(-1i*2*pi* (LineSpacing') * LineSpacing * betas(l)/ (N+1));
    BetaMaps( :, conjLocations ) =  exp(1i*2*pi* (LineSpacing') * LineSpacing * betas(l)/ (N+1));
end
PremultiplicationFactors(L+1,:) = exp(1i*pi*J*N/(N+1));
PostmultiplicationFactors(L+1,:) = exp(1i*pi*N*K/(N+1));

FastRHS = GPUVectorizedCompute2DPolarDFT( gI,  M, L);
c3time = toc( t ); timeit(@()GPUVectorizedCompute2DPolarDFT( gI,  M, L));
fprintf('Execution time on CPU (myFast, more vectorized) = %1.3f secs, it is %1.2fx faster ! \n', ...
    c3time, ctime/c3time );
disp(['Maximum absolute errors element wise  = ',num2str(max(max(abs(LHS-FastRHS))))]);
% LHS = FastRHS;

%%  Using arrayfire mex, capable on 3 platforms, currently designed to use only the CUDA 
CURRENT_PLAT = PLATFORMS.CUDA;
 
switch(CURRENT_PLAT) 
    case PLATFORMS.CUDA % CUDA Path addition
%          addpath('D:\GitHub\RadonTransform\Polar2DTransformMEXArrayFire\x64\Release\');
%          addpath('D:\GitHub\RadonTransform\SphericalPolar3DTransformMEXArrayFire\x64\Release\');
         addpath('..\MEXArrayFireCUDA');
         gI = gpuArray(double(I));                               % Move array to GPU
    case PLATFORMS.CPU % CPU path addition
%         addpath('D:\GitHub\RadonTransform\TwoDPolarTransformMEXArrayFire\x64\CPU_Debug');
%         addpath('D:\GitHub\RadonTransform\TwoDPolarTransformMEXArrayFire\x64\CPU_Release');
%         gI = I;
    case PLATFORMS.OPENCL
end 
 

%% Execution
M = Size+1;                                    % Number of angles
L = (M-2)/4;                                   % Number of levels
if(rem(M-2,4) ~= 0)
    L = ceil (L);                % use + 1 to compute for 45
end

t = tic();
[realImage, ImagImage]= Polar2DTransformMEXArrayFire(gI,M,L,sizeX, uint32( CURRENT_PLAT));

mexOutput = complex(realImage,ImagImage);
g1time =  toc( t );% gputimeit(@() TwoDRadonArrayFire(gI,M,L,sizeX));       % Measure CPU time
if(CURRENT_PLAT == PLATFORMS.CUDA)
    clear mex;  
end  
fprintf('Execution time on GPU (ArrayFire, fully vectorized! on CUDA) = %1.3f secs, it is %1.1fx faster ! \n', ...
    g1time, ctime/g1time ); 
disp(['Maximum absolute errors element wise  = ',num2str(max(max(abs(LHS-mexOutput))))]); 

%% Superfast Vectorized and Precomputations on GPU !!
t = tic(); 
gI = gpuArray(double(I));                               % Move array to GPU
M = Size+1;                                    % Number of angles
L = (M-2)/4;                                   % Number of levels
hasFortyFiveDegrees = 0;
if(rem(M-2,4) ~= 0)
    hasFortyFiveDegrees = 1;
    L = ceil (L);                % use + 1 to compute for 45
end

n=[0:1:N, -N-1:1:-1]';
J = (0:1:N)';
K = (0:1:N)';
Z = zeros(N+1,N+1, 'like',I);
J = J - N/2;

alphas = cosd((1:1:L)*180/M);
betas = sind((1:1:L)*180/M);

LineSpacing = (-N/2:N/2);
% global BetaMaps;
% global E_ns;
% global PremultiplicationFactors;
% global PostmultiplicationFactors;

BetaMaps = zeros ( N+1,2*L*(N+1), 'like',I);
E_ns = zeros ( L,2*N+2, 'like',I);
PremultiplicationFactors = zeros (L+1, N+1,'like',I);
PostmultiplicationFactors = zeros (L+1, N+1,'like',I);

lastIndex = 0;
for l=1:L
    E_ns (l,:) =  exp(-1i*pi*alphas(l)*n.^2/(N+1));                  % Sequence as defined in the paper
    PremultiplicationFactors(l,:) = exp(1i*pi*J *N*alphas(l)/(N+1));
    PostmultiplicationFactors(l,:) = exp(1i*pi*alphas(l)*N*K/(N+1));
    firstLocations =  (1 + lastIndex): lastIndex + (N+1);
    conjLocations = (1 + lastIndex+ (N+1)): lastIndex + 2*(N+1);
    lastIndex = lastIndex + 2*(N+1);
    BetaMaps( :,firstLocations ) =  exp(-1i*2*pi* (LineSpacing') * LineSpacing * betas(l)/ (N+1));
    BetaMaps( :, conjLocations ) =  exp(1i*2*pi* (LineSpacing') * LineSpacing * betas(l)/ (N+1));
end
PremultiplicationFactors(L+1,:) = exp(1i*pi*J*N/(N+1));
PostmultiplicationFactors(L+1,:) = exp(1i*pi*N*K/(N+1));


FastRHSGPU = GPUVectorizedCompute2DPolarDFT( gI,  M, L);
gtime = toc( t ); gputimeit(@()GPUVectorizedCompute2DPolarDFT( gI,  M, L));       % Measure CPU time
% fprintf('Execution time on GPU (mySuperFast, more vectorized) = %1.3f secs, it is %1.1fx faster ! \n', ...
%     gtime, ctime/gtime );
fprintf('Execution time on GPU (mySuperFast, more vectorized) = %1.3f secs, it is %1.1fx faster ! \n', ...
    gtime, ctime/gtime );
disp(['Maximum absolute errors element wise  = ',num2str(max(max(abs(LHS-FastRHSGPU))))]); 
