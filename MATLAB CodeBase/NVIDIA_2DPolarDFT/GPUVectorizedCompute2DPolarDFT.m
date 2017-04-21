function [ PolarGrid ] = GPUVectorizedCompute2DPolarDFT( inputImage,  noOfAngles,noOfLevels)
%GPUVectorizedCompute2DPolarDFT does a super fast computations
%   for a radon transform
I = inputImage;
global BetaMaps;
global E_ns;
global PremultiplicationFactors;
global PostmultiplicationFactors;
global Z;

[sizeX,~] =  size(I);
N = sizeX -1;      % N is always even
M = noOfAngles;    % M is also even
L = noOfLevels;    % This precomputed and given here too

PolarGrid = zeros(M, N+1, 'like',I);     %  Polar grid: No of angles vs. Radial data

%% Actual loop
lastIndex = 0;
for l=1:1:L  % For each level
    angle = l*180/M;
    
    gpuEn  = (E_ns (l,:).');  %------------- GPU
    gpuPreMulFactor = (PremultiplicationFactors(l,:).');
    gpuPostMulFactor = (PostmultiplicationFactors(l,:).');
    
    firstLocations =  (1 + lastIndex): lastIndex + (N+1);
    conjLocations = (1 + lastIndex+ (N+1)): lastIndex + 2*(N+1);
    lastIndex = lastIndex + 2*(N+1);
    BetaMap =  BetaMaps(  :, firstLocations) ;
    BetaMapConj = BetaMaps(  :,conjLocations );
    
    gpu_0_90PreMul = PremultiplicationFactors(L+1,:);
    gpu_0_90PostMul = PostmultiplicationFactors(L+1,:);
    
    % X axis Scaling or row wise scaling
    F_x_alpha_l = GPUVectorizedFrFT_Centered( I.',Z, gpuEn ,gpuPreMulFactor ,gpuPostMulFactor);
    F_x_alpha_l = F_x_alpha_l.';
     
    if (l == 1) % For first pass gather lines at 0 and 90
        NintyLine = F_x_alpha_l(:, N/2+1);          % Getting the column of 90 degrees angle
        col = bsxfun(@times,gpu_0_90PostMul.',fft(bsxfun(@times,NintyLine,gpu_0_90PreMul.')));
        PolarGrid (M/2+1,:)   =  ( col);
    end
    
    PolarGrid (1+l, :) =  sum ( bsxfun(@times,(F_x_alpha_l),BetaMap));
    PolarGrid (M+1-l, :) = fliplr( sum ( bsxfun(@times,(F_x_alpha_l),BetaMapConj)));
      
    if (angle == 45)   % If angle is 45 we have already computed it above
        continue;
    end 
     
    % Y axis Scaling          
    F_y_alpha_l = GPUVectorizedFrFT_Centered((I),Z, gpuEn ,gpuPreMulFactor ,gpuPostMulFactor);
    
    if (l == 1)         % For first pass gather lines at 0 and 90
        ZeroLine = F_y_alpha_l( N/2+1,:);               % Getting the row of 0 degrees angle
        PolarGrid (1,:) =  bsxfun( @times,gpu_0_90PostMul,fft(bsxfun(@times,ZeroLine,gpu_0_90PreMul)) ); % At 0 degrees
    end
    
    PolarGrid (M/2+1+l, :) = sum ( bsxfun(@times,F_y_alpha_l.',BetaMapConj) );
    PolarGrid (M/2+1-l, :) = sum ( bsxfun(@times,F_y_alpha_l.',BetaMap) );
    
end 
end 

