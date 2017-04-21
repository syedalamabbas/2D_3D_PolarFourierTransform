

clear all
close all
clc

Sizes =  15:2: 799;

ctimeDirect = 10^12*zeros(1,length(Sizes));
ctimeNonVector = zeros(1,length(Sizes));
ctimeNonHalfVectorized = zeros(1,length(Sizes));
ctimeCPUVectorized = zeros(1,length(Sizes));
ctimeGPUVectorized = zeros(1,length(Sizes));
global BetaMaps; 
global E_ns;
global PremultiplicationFactors;
global PostmultiplicationFactors;
global Z ;

for k =1:length(Sizes)
    disp(['Iteration-', num2str(k)]);
    Size = Sizes(k);
    N = Size -1;
    M = Size +1;
    I = rand(Size);
    if(Size <= 160 )
        ctimeDirect(k) = timeit(@()DirectBruteForce2D_DFT( I,  M));
    end
    ctimeNonVector(k) = timeit(@()Compute2DPolarDFT( I,  M));
    ctimeNonHalfVectorized(k) =   timeit(@()VectorizedCompute2DPolarDFT( I,  M));
    %%
    L = (M-2)/4;                                   % Number of levels
    if(rem(M-2,4) ~= 0)
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
    
    ctimeCPUVectorized(k) =  timeit(@()GPUVectorizedCompute2DPolarDFT( I,  M, L));
    
    %% gpuArray(I)
    I = gpuArray(I);
     L = (M-2)/4;                                   % Number of levels
    if(rem(M-2,4) ~= 0)
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
    
    ctimeGPUVectorized(k) = gputimeit(@()GPUVectorizedCompute2DPolarDFT( I, M, L));
end

%% Plot the timining of the runs
N = 393;

% Some averaging to make the graph look pretty, this is because we are not
% cleaning the GPU variables properly , look into this ?
Array1 = ctimeGPUVectorized;
Array2 = reshape(Array1,3,[]);
Array3  = sum(Array2,1)./ size(Array2,1);
ctimeGPUVectorized = reshape(repmat(Array3, size(Array2,1),1 ) , 1,[]);

figure,
semilogy(Sizes(1:N),ctimeDirect(1:N), 'r', 'LineWidth', 2.3)
hold on
semilogy(Sizes(1:N),ctimeNonVector(1:N),'b','LineWidth', 2.3)
semilogy(Sizes(1:N),ctimeNonHalfVectorized(1:N),'g','LineWidth', 2.3)
semilogy(Sizes(1:N),ctimeCPUVectorized(1:N),'m','LineWidth', 2.3)
semilogy(Sizes(1:N),ctimeGPUVectorized(1:N),'k','LineWidth', 2.3)
hold off
xlabel('(N+1), is the number of samples for (N+1)^2  data')
ylabel('Average running time (seconds)')
legend('Direct Brute Force', 'Fast+ non-vectorized', 'Fast+ vectorized' , 'Fast + vectorized + precomputation', 'GPU + vectorized + precomputation')
grid on
axis tight
print -dpdf ImagesRunningTime 