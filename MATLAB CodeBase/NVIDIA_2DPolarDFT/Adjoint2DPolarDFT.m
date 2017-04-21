function [ ImageAdjoint ] = Adjoint2DPolarDFT( Polar_Grid )
%ADJOINT2DPOLARDFT computes the adjoint linear operation as described in
%the paper by Syed Alam Abbas

%====================================================================
% This function performs a Polar ADJOINT transform on a 2D signal X
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
%     N=16;                                            % Always even
%     M = N +2 ;                                       % Always even
%     X=randn(N+1, N+1) + 1i*randn(N+1, N+1);          % Complex matrix for Polar Grid
%     Y=randn(M, N+1) + 1i*randn(M, N+1);              % Complex matrix of Polar Grid
%     Ax=Compute2DPolarDFT( X,  M );                   % Output M x N+1    polar grid
%     AHy=Adjoint2DPolarDFT(Y);                        % Output  N+1 x N+1   adjoint image
%     prod1 = dot(Ax,Y);   
%     prod2 = dot(X,AHy.');                            % Dont know why we need an additional element wise transpose, also used in inverse operation 
%     fprintf('\n The result of  (<Ax,y> - <x, A^H y>) comparison is :')
%     disp(abs(sum(sum(prod1)) -(sum(sum(prod2)))));
%    Output Like ----------->>>>  The result of  (<Ax,y> - <x, A^H y>) comparison is :   3.6265e-11
% Written on June 29th, 2015 by Syed Alam Abbas.
%====================================================================

[M, sizeN_plus1 ] = size(Polar_Grid);
N = sizeN_plus1 - 1;

L = (M-2)/4;  % Number of levels
if(rem(M-2,4) ~= 0)
    L = ceil (L);                   % use + 1
end

lineSpacing = -N/2:N/2;
gridSpacing = lineSpacing.' * lineSpacing;

ImageAdjoint = zeros (sizeN_plus1, sizeN_plus1); % Suppose to be Complex Adjoint image

for l=1:1:L  % For each level
    
    angle = l*180/M;
    alpha_l = cosd(angle);
    beta_l  = sind(angle);
    
    beta_factor = exp(2*pi*1i* beta_l * gridSpacing  / sizeN_plus1); % observe the scale change, it is +ve
    
    % Reversing the forward Polar Grid operation step by step
    
    %% X-axis FrFT Grid
    
    % computing variable scale FrFT  first, gather the  lines oriented at angle from x-axis
    line1 = Polar_Grid (1+l, :);
    line2 = fliplr(Polar_Grid (M+1-l, :));
    
    FrFTVarBlock = getBlockData2ComplementaryLines( line1, line2,beta_factor); % they have the same scale so multiplication FrFT alpha is the same
    
    % computing uniform scale FrFT second
    FrFTUniformBlock_x = VectorizedFrFT_Centered(FrFTVarBlock.' , -alpha_l);       % observe the scale change , X axis scaling in reverse
    
    ImageAdjoint = ImageAdjoint + FrFTUniformBlock_x;                              % collecting contribution from every level
    
    
    if (angle == 45)
        continue;
    end
    
    %% Y-axis FrFT Grid
    line1 = Polar_Grid (M/2+1-l,:);
    line2 = Polar_Grid (M/2 + 1+ l, :);
    FrFTVarBlock = getBlockData2ComplementaryLines( line1, line2,beta_factor); % They have the same scale , so the multiplication of alpha values is the same
    
    % computing uniform scale FrFT second
    FrFTUniformBlock_y = VectorizedFrFT_Centered(FrFTVarBlock.' , -alpha_l);     % observe the scale change, Y axis scaling in reverse
    
    ImageAdjoint = ImageAdjoint + FrFTUniformBlock_y.';                         % collecting contribution from every level
    
    %% computing for zero and ninety degrees seperately
    
    if (l == 1)
        ZeroLine = Polar_Grid (1,:);
        NinetyLine = Polar_Grid(M/2+1,:);
        
        FrFTVarBlock1 =  ones(sizeN_plus1,1)* NinetyLine;        % beta factors is all ones
        FrFTVarBlock2 =  ones(sizeN_plus1,1) * ZeroLine;    % beta factors is all ones
        
        % computing uniform scale FrFT second
        FrFTUniformBlock_x = VectorizedFrFT_Centered(FrFTVarBlock2.' , -1);          % Zero
        FrFTUniformBlock_y = VectorizedFrFT_Centered(FrFTVarBlock1.' , -1);
        
        ImageAdjoint = ImageAdjoint + FrFTUniformBlock_x + FrFTUniformBlock_y.';   % collecting contribution from every level
    end
    
end

    function BlockDataFrom2Lines = getBlockData2ComplementaryLines(lineData1, lineData2, beta_factor)
        tiledLineData =   ones ( sizeN_plus1, 1) * lineData1;
        tiledLineDataConj =   ones ( sizeN_plus1, 1) * lineData2;
        
        FrFTVarBlock1 = beta_factor .* tiledLineData;
        FrFTVarBlock2 = conj(beta_factor) .* tiledLineDataConj;
        
        BlockDataFrom2Lines = FrFTVarBlock1 + FrFTVarBlock2; % they have the same scale so multiplication FrFT alpha is the same
        
    end
end

