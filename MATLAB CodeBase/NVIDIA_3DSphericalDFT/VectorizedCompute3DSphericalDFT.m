function [ SphericalGrid ] = VectorizedCompute3DSphericalDFT( inputImage,  noOfAnglesTheta, noOfAnglesPhi ,Debug3D_DFT )
%COMPUTE2DRADONTRANSFORM Transforms the fucntion using the prescribed
% scheme as discussed in the paper by Syed Alam Abbas
%====================================================================
% Written on October 7th, 2015 by Syed Alam Abbas.
% Updated on 12/9/2015 , fixed a X-Y switch bug
%====================================================================

%% Initializing variables some are declared global to be used by internal functions
I = inputImage;
global N;
global d_NoOfElements;
global lineSpacing;


[sizeX,~,~] =  size(I);
d_NoOfElements = sizeX;
N = sizeX -1;      % N is always even
lineSpacing = -N / 2: N / 2;
M = noOfAnglesPhi;    % M is also even
K = noOfAnglesTheta;

%% Computing number of levels
% No of levels theta                           % These computed on fly
P = (K-2)/4;                                   % Number of levels
if(rem(K-2,4) ~= 0)
    P = ceil (P);                           % use + 1 to compute for 45 but skip for other dimension
end

% No of levels Phi
Q = (M-2)/4;                                   % Number of levels
if(rem(M-2,4) ~= 0)
    Q = ceil (Q);                           % use + 1 to compute for 45 but skip for other dimension
end

%% These are the required additional variables
SphericalGrid = zeros(K, M, N+1, 'like',I);                % Angle phi vs  Polar slices: No of angles vs. Radial data
% ReorderedImage_OperateColumns;                 % Reordered 3D Image at the start of the FrFT operation
% FrFT1D_Uniform_Image3D;                        % 3D Image obtained after passing 1st stage of FrFT, all columns have now been operated on
% Slice2D, Slice2D_Conj;                         % 2D Slices obtained at the 2nd stage of FrFT
% Line1, Line1_Conj, Line2, Line2_Conj;		     % 1D Lines at the end of the 3rd stage of FrFT, which are final
% alpha_factor, beta_factor, gamma_factor;       % Scaling factors in X-axis , Y-axis and Z-axis respectively, which change depending on the computation of the block

debug = 0;                                       % Close this flag to stop debugging, requires precomputed DirectBruteForce3D_DFT,set this flag to 0 and send null when not debugging

%% Processing
for q=1:Q  % For each level of Polar slices
    anglePhi = q*180/M;
    
    %% X oriented pair of Polar slices, (1+q, K+1-q)
    for p = 1:P
        angleTheta = p*180/K;
        
        % XX block  -- Concentric rectangles in YZ tiled along X -axis  (1+p, M+1-p)
        alpha_factor = cosd(angleTheta) * cosd(anglePhi);											         % Scaling needed as defined in the paper
        beta_factor = cosd(angleTheta) * sind(anglePhi);
        gamma_factor = sind(angleTheta);
        
        ReorderedImage_OperateColumns = I;                                  %  No Swap, X-axis as columns now
        FrFT1D_Uniform_Image3D = VectorizedFrFT_Centered(ReorderedImage_OperateColumns,alpha_factor);  % X-axis uniform scaling
        
        if (q == 1 && p == 1)         %% Computing Polar Slice at anglephi = 90, only needed to be computed once
            Central2D_YZSlice2D = FrFT1D_Uniform_Image3D(N/2+1, :, :);
            Central2D_YZSlice2D = squeeze(Central2D_YZSlice2D);
            properOrientedSlice = (Central2D_YZSlice2D.');          %% Verified match, yes
            Polar2D = VectorizedCompute2DPolarDFT( properOrientedSlice,  K );
            SphericalGrid(:, M/2+1, :) = Polar2D;
            if(debug)
                disp('This is Polar Slice at angle phi = 90');
                squeeze(Debug3D_DFT(:, M/2+1, :) - SphericalGrid(:, M/2+1, :))
            end
        end
        
        FrFT1D_Uniform_Image3D = permute (FrFT1D_Uniform_Image3D, [ 2 1 3]);
        % Y axis differential scaling
        [Slice2D, Slice2D_Conj] = ComputeNonVector_NonUniformFrFT_ComplementaryPlanes(FrFT1D_Uniform_Image3D, beta_factor);
        
        % Z axis differential scaling
        [Line1, Line1_Conj]= ComputeFinal_2ComplementaryLines(Slice2D.', gamma_factor);
        [Line2, Line2_Conj]= ComputeFinal_2ComplementaryLines(Slice2D_Conj.',gamma_factor);
        
        %% Finally collecting lines for XX block
        SphericalGrid(p+1, q+1, :) = Line1;
        SphericalGrid(p+1, M+1 - q, :) = conj(Line2_Conj);
        SphericalGrid(K+1 - p, q+1, :) = conj(Line1_Conj);											%% 1 More Swap to match !
        SphericalGrid(K+1 - p, M+1 - q, :) = Line2;
        %% Used only for debugging
        if(debug)
            disp('This is XX block solution');
            squeeze(Debug3D_DFT(p+1, q+1, :) - SphericalGrid(p+1, q+1, :))
            squeeze(Debug3D_DFT(p+1, M+1 - q, :)-SphericalGrid(p+1, M+1 - q, :))
            squeeze(Debug3D_DFT(K+1 - p, q+1, :)-SphericalGrid(K+1 - p, q+1, :) )
            squeeze(Debug3D_DFT(K+1 - p, M+1 - q, :)-SphericalGrid(K+1 - p, M+1 - q, :) )
        end
        
        %% If theta = 45, skip the next block which is redundant
        if (angleTheta == 45 )
            continue;
        end
        
        %% XZ block  -- Concentric rectangles in
        alpha_factor = sind(angleTheta) * cosd(anglePhi);											         % Scaling needed as defined in the paper
        beta_factor  = sind(angleTheta) * sind(anglePhi);
        gamma_factor = cosd(angleTheta);
        
        ReorderedImage_OperateColumns = permute (I , [3 2 1]);                            %  1st Swap, Z-axis as columns now, operate along Z dimension
        FrFT1D_Uniform_Image3D = VectorizedFrFT_Centered(ReorderedImage_OperateColumns,gamma_factor);
        
        if (q == 1 && p == 1)         %% Computing One special Polar Slice at angletheta = 0 , only need to be computed once
            Central2D_XYSlice2D = FrFT1D_Uniform_Image3D(N / 2+1, :, :);
            Central2D_XYSlice2D = squeeze(Central2D_XYSlice2D);
            properOrientedSlice = (Central2D_XYSlice2D);          %% Verified match  , yes
            Polar2D = VectorizedCompute2DPolarDFT(properOrientedSlice, M);
            SphericalGrid(1, :, :) = Polar2D;
            if(debug)
                disp('This is Polar Slice at angle theta = 0');
                squeeze(Debug3D_DFT(1, :, :) - SphericalGrid(1, :, :))
            end
        end
        FrFT1D_Uniform_Image3D = permute (FrFT1D_Uniform_Image3D, [ 2 1 3]);
        
        % Y axis differential scaling
        [Slice2D, Slice2D_Conj] = ComputeNonVector_NonUniformFrFT_ComplementaryPlanes(FrFT1D_Uniform_Image3D, beta_factor);
        
        % X axis differential scaling
        [Line1, Line1_Conj]= ComputeFinal_2ComplementaryLines(Slice2D.', alpha_factor);
        [Line2, Line2_Conj]= ComputeFinal_2ComplementaryLines(Slice2D_Conj.',alpha_factor);
        %             if(p == 1 && q == 1 )
        %                 if(p == 2 && q == 1)
        %         if(p == 2 && q == 2)
        %         if(p == 1 && q == 2)
        %         if((p == 2 || p == 1) && (q == 1 || q == 2))
        SphericalGrid(K/2+1 - p, 1+q, :)      = Line1;
        SphericalGrid(K/2+1 - p, M+1 - q, :)  = (Line1_Conj);
        SphericalGrid(K/2+1 + p, q+1, :)      = (Line2_Conj);
        SphericalGrid(K/2+1 + p, M+1 - q, :)  = (Line2);
        %         end
        
        %% Used only for debugging
        if(debug)
            disp('This is XZ block solution');
            squeeze(Debug3D_DFT(K/2+1 - p, 1+q, :) - SphericalGrid(K/2+1 - p, 1+q, :))
            squeeze(Debug3D_DFT(K/2+1 - p, M+1 - q, :)-SphericalGrid(K/2+1 - p, M+1 - q, :))
            squeeze(Debug3D_DFT(K/2 + p +1, q+1, :)-SphericalGrid(K/2 + p +1, q+1, :) )
            squeeze(Debug3D_DFT(K/2 + p +1, M+1 - q, :)-SphericalGrid(K/2 + p +1, M+1 - q, :) )
        end
        
        %% If Phi = 45, skip the next block which is redundant
        if (anglePhi == 45)   % If angle is 45 we have already computed it above
            continue;
        end
        
        %% YZ block -- Concentric rectangles in
        alpha_factor = sind(angleTheta) * sind(anglePhi);		   % Scaling needed as defined in the paper
        beta_factor  = sind(angleTheta) * cosd(anglePhi);
        gamma_factor = cosd(angleTheta);                           % see it is the same scale so use the previously computed grid FrFT 1D
        
        % Y axis differential scaling
        [Slice2D, Slice2D_Conj] = ComputeNonVector_NonUniformFrFT_ComplementaryPlanes(FrFT1D_Uniform_Image3D, beta_factor);
        
        % X axis differential scaling
        [Line1, Line1_Conj]= ComputeFinal_2ComplementaryLines(Slice2D.', alpha_factor);
        [Line2, Line2_Conj]= ComputeFinal_2ComplementaryLines(Slice2D_Conj.',alpha_factor);
        
        SphericalGrid(K / 2+1 - p, M / 2+1 - q, :) = Line1;
        SphericalGrid(K / 2+1 - p, M / 2+1 + q, :) = Line1_Conj;
        SphericalGrid(K/ 2+1 + p, M/ 2+1 - q, :) = Line2_Conj;
        SphericalGrid(K/ 2+1 + p, M / 2+1 + q, :) = Line2;
        
        %% Used only for debugging
        if(debug)
            disp('This is YZ block solution');
            squeeze(Debug3D_DFT(K / 2+1 - p, M / 2+1 - q, :) - SphericalGrid(K / 2+1 - p, M / 2+1 - q, :))
            squeeze(Debug3D_DFT(K / 2+1 - p, M / 2+1 + q, :)-SphericalGrid(K / 2+1 - p, M / 2+1 + q, :))
            squeeze(Debug3D_DFT(K/ 2+1 + p, M/ 2+1 - q, :)-SphericalGrid(K/ 2+1 + p, M/ 2+1 - q, :) )
            squeeze(Debug3D_DFT(K/ 2+1 + p, M / 2+1 + q, :)-SphericalGrid(K/ 2+1 + p, M / 2+1 + q, :) )
        end
    end
    
    %% If Phi = 45, skip the next block which is redundant
    if (anglePhi == 45)   % If angle is 45 we have already computed it above
        continue;
    end
    
    %%  Y oriented pair of Polar slices, (d_NoOfAnglesPhi / 2 - q, d_NoOfAnglesPhi / 2 + q)
    for p = 1:P
        angleTheta = p*180/K;
        %% YY block
        alpha_factor = cosd(angleTheta) * sind(anglePhi);				 %% Scaling needed as defined in the paper
        beta_factor = cosd(angleTheta) * cosd(anglePhi);
        gamma_factor = sind(angleTheta);
        
        ReorderedImage_OperateColumns =   permute (I , [2 1 3]);
        
        % Y-axis uniform scaling
        FrFT1D_Uniform_Image3D = VectorizedFrFT_Centered(ReorderedImage_OperateColumns,beta_factor);
        
        if (q == 1 && p == 1)         %% Computing One special Polar Slice at  anglephi = 0 , only need to be computed once
            Central2D_YZSlice2D = FrFT1D_Uniform_Image3D(N / 2+1, :, :);
            Central2D_YZSlice2D = squeeze(Central2D_YZSlice2D);
            properOrientedSlice = (Central2D_YZSlice2D.');          %% Verified match  yes
            Polar2D = VectorizedCompute2DPolarDFT( properOrientedSlice, K );
            SphericalGrid(:, 1, :) = Polar2D;
            if(debug)
                disp('This is Polar Slice at angle phi = 0');
                squeeze(Debug3D_DFT(:, 1, :) - SphericalGrid(:, 1, :))
            end
            
            % Special operation
            SpecialLineZ = Polar2D(K / 2+1, :);
            TiledZ = repmat(SpecialLineZ,M,1);
            SphericalGrid(K/2+1, :, :) = TiledZ;     %% This line is common to all grids !!! VERY VERY Special
            if(debug)
                disp('This VERY VERY Special line (common to so many polar girds) at angle phi = 0 and angle theta = 90');
                squeeze(Debug3D_DFT(K/2+1, :, :) - SphericalGrid(K/2+1, :, :))
            end
        end
        FrFT1D_Uniform_Image3D = permute (FrFT1D_Uniform_Image3D, [ 2 1 3]);
        %% X axis differential scaling
        [Slice2D, Slice2D_Conj] = ComputeNonVector_NonUniformFrFT_ComplementaryPlanes(FrFT1D_Uniform_Image3D, alpha_factor);
        
        %% Z axis differential scaling
        [Line1, Line1_Conj]= ComputeFinal_2ComplementaryLines(Slice2D.',gamma_factor );
        [Line2, Line2_Conj]= ComputeFinal_2ComplementaryLines(Slice2D_Conj.',gamma_factor);
        
        SphericalGrid(p+1, M/ 2+1 - q, :) = Line1;
        SphericalGrid(p+1, M/ 2+1 + q, :) = Line2;
        SphericalGrid(K+1 - p, M/ 2+1 - q, :) = conj(Line1_Conj);
        SphericalGrid(K+1 - p, M/ 2+1 + q, :) = conj(Line2_Conj);
        
        %% Used only for debugging
        if(debug)
            disp('This is YY block solution');
            squeeze(Debug3D_DFT(p+1, M/ 2+1 - q, :) - SphericalGrid(p+1, M/ 2+1 - q, :))
            squeeze(Debug3D_DFT(p+1, M/ 2+1 + q, :)-SphericalGrid(p+1, M/ 2+1 + q, :))
            squeeze(Debug3D_DFT(K+1 - p, M/ 2+1 - q, :)-SphericalGrid(K+1 - p, M/ 2+1 - q, :) )
            squeeze(Debug3D_DFT(K+1 - p, M/ 2+1 + q, :)-SphericalGrid(K+1 - p, M/ 2+1 + q, :) )
        end
    end
end
    function [Slice2D, Slice2D_Conj] = ComputeNonVector_NonUniformFrFT_ComplementaryPlanes(FrFT1D_Image_3D, ColumnScale)
        lineSpacing_Cube_TiledLevel = repmat(lineSpacing'*lineSpacing,[1 1 d_NoOfElements]);
        Column_ScaleFactor = exp (-2*1i * pi * lineSpacing_Cube_TiledLevel * ColumnScale / d_NoOfElements);
        
        %% Consider multiplication of complex numbers A = (a+ib); B = (c+id) four double real multiplications
        realRealPart = real(FrFT1D_Image_3D) .* real(Column_ScaleFactor);             %% ac
        realImagPart = real(FrFT1D_Image_3D) .* imag(Column_ScaleFactor);             %% ad
        imagRealPart = imag(FrFT1D_Image_3D) .* real(Column_ScaleFactor);             %% bc
        imagImagPart = imag(FrFT1D_Image_3D) .* imag(Column_ScaleFactor);             %% bd
        
        Slice2D = sum((realRealPart - imagImagPart)+ 1i* (realImagPart + imagRealPart));       %% sum(flip(FrFT_Image_X_Cube, 0) *BetaFactor);
        Slice2D_Conj = sum((realRealPart + imagImagPart) + 1i*( imagRealPart - realImagPart));  %%  sum(flip(FrFT_Image_X_Cube, 0) * conjg(BetaFactor));
        
        Slice2D = squeeze(Slice2D);                                                 %% Making it 2D without changing the order of data
        Slice2D_Conj =squeeze(Slice2D_Conj);
    end

    function [Line1D, Line1D_Conj] = ComputeFinal_2ComplementaryLines(FrFTFrFT_Image2D, ColumnScale)
        lineSpacing_Square_TiledLevel = lineSpacing'*lineSpacing;
        Column_ScaleFactor = exp(-2*1i * pi * lineSpacing_Square_TiledLevel * ColumnScale / d_NoOfElements);
        
        %% Consider multiplication of complex numbers A = (a+ib); B = (c+id) four double real multiplications
        realRealPart = real(FrFTFrFT_Image2D) .* real(Column_ScaleFactor);             %% ac
        realImagPart = real(FrFTFrFT_Image2D) .* imag(Column_ScaleFactor);             %% ad
        imagRealPart = imag(FrFTFrFT_Image2D) .* real(Column_ScaleFactor);             %% bc
        imagImagPart = imag(FrFTFrFT_Image2D) .* imag(Column_ScaleFactor);             %% bd
        
        Line1D = sum((realRealPart - imagImagPart)+ 1i* (realImagPart + imagRealPart));
        Line1D_Conj = sum((realRealPart + imagImagPart) + 1i*( imagRealPart - realImagPart));
    end
end