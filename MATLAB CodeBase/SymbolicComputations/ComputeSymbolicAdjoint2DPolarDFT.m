function [ ImageAdjoint ] = ComputeSymbolicAdjoint2DPolarDFT( Polar_Grid, alphas, betas, one_alpha )
%ADJOINT2DPOLARDFT computes the adjoint linear operation as described in
%the paper

%====================================================================
% This function performs a Polar ADJOINT transform symbolically  on a 2D signal X 
% given on the polar coordinate system. If X is M x (N+1), the output will have (N+1) x (N+1) values. 
% The algorithm applied here uses the Fast Fractional Fourier Transform. 
%   
% Synopsis: Y=Adjoint2DPolarDFT(X)
%
% Inputs -  X     M x (N+1) matrix in polar grid, (N is assumed to be even)           
% Outputs - Y      (N+1) x (N+1) matrix (full Cartesian grid)
%
% Example: 
% 
%   The following is a way to verify that this works as an adjoint -
%   choosing random X and Y, one must obtain that 
% The adjoint of a transformation is a unique transformation T* so that 
% <Tx,y> = <x, T*y> for every x and y
% <y,Ax> = <x,adj(A)y>'
%  
%   N=16;                                            % Always even
%   M = N +2 ;                                       % Always even
%   X=randn(N+1, N+1)+sqrt(-1)*randn(N+1, N+1);      % Complex matrix of Polar Grid
%   Y=randn(M, N+1)+sqrt(-1)*randn(M, N+1);          % Complex matrix
%   AX=Compute2DPolarDFT( X,  M );                   % output M x N+1       polar   grid
%   AtY=Adjoint2DPolarDFT(Y);                        % Output  N+1 x N+1   adjoint image
%   disp(abs( sum(sum(Y'.*AX)) -conj(sum(sum(X'.*AtY)))));       
% 
% Written on June 29th, 2015 by Syed Alam Abbas.
%====================================================================


[M, sizeN_plus1 ] = size(Polar_Grid);
N = sizeN_plus1 - 1;

L = (M-2)/4;  % Number of levels
if(rem(M-2,4) ~= 0)
    L = ceil (L);                   % use + 1
end

spacing = -N/2:N/2;
gridSpacing = spacing.' * spacing;

ImageAdjoint = zeros (sizeN_plus1, sizeN_plus1); % Complex Adjoint image

for l=1:1:L  % For each level

    angle = l*180/M;
    alpha_l = alphas(l); cosd(angle);
    beta_l  = betas(l); sind(angle);
    
    beta_factor = exp(2*pi*1i* beta_l* gridSpacing  / sizeN_plus1); % observe the scale change, it is +ve
    
    % Reversing the forward Polar Grid operation step by step
    
    %% X-axis FrFT Grid
    
%     % computing variable scale FrFT  first
    lineData = Polar_Grid (1+l, :);   % gather the  line oriented at angle from x-axis
    lineData2 = Polar_Grid (M+1-l, :);   % gather the  line oriented at angle from x-axis but in opposite direction
    
    tiledLineData =  ones ( sizeN_plus1, 1) * lineData;      
    tiledLineDataConj =   ones ( sizeN_plus1, 1) * fliplr(lineData2);
    
    FrFTVarBlock1 = beta_factor .* tiledLineData;          % Switched
    FrFTVarBlock2 = conj(beta_factor) .* tiledLineDataConj;
    
    FrFTVarBlock = FrFTVarBlock1 + FrFTVarBlock2; % they have the same scale so multiplication FrFT alpha is the same
    
    % computing uniform scale FrFT second
    % X axis Scaling
    for p = 1: sizeN_plus1      % For row
        for q = 1: sizeN_plus1   % For columns
            FrFTUniformBlock_x(p,q) = SymbolicFrFT_Uniform( FrFTVarBlock(p,:), -alpha_l, spacing(q), sizeN_plus1 );   % observe the scale change , X axis scaling in reverse
        end
    end
    ImageAdjoint = ImageAdjoint + FrFTUniformBlock_x.';                            % collecting contribution from every level
    
    if (angle == 45)
        continue;
    end
    
    %% Y-axis FrFT Grid
    lineDataY = Polar_Grid (M/2+1-l,:);
    lineDataY2 = Polar_Grid (M/2 + 1+ l, :);
    
    tiledLineDataY =  lineDataY' * ones(1, sizeN_plus1) ;
    tiledLineDataConjY =  lineDataY2' * ones(1, sizeN_plus1)  ;
    
    FrFTVarBlock1 = beta_factor .* tiledLineDataY;         
    FrFTVarBlock2 = conj(beta_factor) .* tiledLineDataConjY;       
    
    FrFTVarBlock = FrFTVarBlock1 + FrFTVarBlock2; % They have the same scale , so the multiplication of alpha values is the same
    
    % computing uniform scale FrFT second 
    
     % Y axis Scaling
    for p = 1: sizeN_plus1      % For row
        for q = 1: sizeN_plus1   % For columns
            FrFTUniformBlock_y(q,p) = SymbolicFrFT_Uniform( FrFTVarBlock(:,p).', -alpha_l, spacing(q), sizeN_plus1 );      
        end
    end

    ImageAdjoint = ImageAdjoint + FrFTUniformBlock_y.';                           % collecting contribution from every level
    
    %% computing for zero and ninety degrees seperately
    
    if (l == 1)
        zeroLine = Polar_Grid (1,:);
         NinetyLine = Polar_Grid(M/2+1,:);
        
%         tiledLineDataOrtho = zeroLine.' * ones(1, sizeN_plus1)  ;
%         tiledLineDataConjOrtho = ones(sizeN_plus1, 1) * NinetyLine;
        
        tiledLineDataOrtho = NinetyLine.' * ones(1, sizeN_plus1)  ;
        tiledLineDataConjOrtho = ones(sizeN_plus1, 1) * zeroLine;
        
        
        FrFTVarBlock1 =  tiledLineDataOrtho;
        FrFTVarBlock2 = tiledLineDataConjOrtho;    % beta factors is all ones ???
        
        
        
        % X axis Scaling
        for p = 1: sizeN_plus1      % For row
            for q = 1: sizeN_plus1   % For columns
                FrFTUniformBlock_x(p,q) = SymbolicFrFT_Uniform( FrFTVarBlock2(p,:), -one_alpha, spacing(q), sizeN_plus1 );
            end
        end
        
        % Y axis Scaling
        for p = 1: sizeN_plus1      % For row
            for q = 1: sizeN_plus1   % For columns
                FrFTUniformBlock_y(q,p) = SymbolicFrFT_Uniform( FrFTVarBlock1(:,p).', -one_alpha, spacing(q), sizeN_plus1 );
            end
        end

        ImageAdjoint = ImageAdjoint + FrFTUniformBlock_x.' + FrFTUniformBlock_y.';                   % collecting contribution from every level
        
    end
    
end
end

