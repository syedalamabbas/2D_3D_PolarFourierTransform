clear ;
clc;

for i = 1: 1
    disp(['iteration', num2str(i)])
    %% Sizes
    N = 3024;
    M = 1;
    a = complex(randn(N,M ),randn(N,M ));   % Data input
    b = complex(randn(N,M ),randn(N,M ));randn(N,1);                                % Filter input
    
    %% CPU Non-Vectorized
    c = fastConvolution(a,b);                      % Calculate output
    ctime = timeit(@()fastConvolution(a,b));       % Measure CPU time
    disp(['Execution time on CPU(non-vector) (in seconds)= ',num2str(ctime)]);
    
    
    %% Special 
     
    %% GPU Non-vectorized
    ga = gpuArray(a);                              % Move array to GPU
    gb = gpuArray(b);                              % Move filter to GPU
    gc = fastConvolution(ga,gb);                   % Calculate on GPU
    gtime = gputimeit(@()fastConvolution(ga,gb));  % Measure GPU time
    gerr = max(max(abs(gather(gc)-c)));            % Calculate error
    disp(['Execution time on GPU(non-vector)(in seconds) = ',num2str(gtime)]);
    disp(['Maximum absolute error = ',num2str(gerr)]);
    
    %% CPU Vectorized
    
    c = VectorizedfastConvolution(a,b);                    % Calculate output
    ctime = timeit(@()VectorizedfastConvolution(a,b));     % Measure CPU time
    disp(['Execution time on CPU (Vectorized) (in seconds)= ',num2str(ctime)]);
    
    
    %% GPU Vectorized
    ga = gpuArray(a);                               % Move data to GPU
    gb = gpuArray(b);                               % Move filter to GPU
    gc = VectorizedfastConvolution(ga, gb);                % Calculate on GPU
    gtime = gputimeit(@()VectorizedfastConvolution(ga,gb));% Measure GPU time
    gerr = max(max(abs(gather(gc)-c)));             % Calculate error
    disp(['Execution time on GPU (Vectorized) (in seconds) = ',num2str(gtime)]);
    disp(['Maximum absolute error = ',num2str(gerr)]);
end


