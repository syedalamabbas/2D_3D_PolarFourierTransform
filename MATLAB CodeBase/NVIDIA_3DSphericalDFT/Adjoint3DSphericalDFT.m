function [ Image3DAdjoint ] = Adjoint3DSphericalDFT( SphericalGrid , DebugAdjoint_3DDFT)
%ADJOINT3DSPHERICALDFT computes the adjoint linear operation for 3D data as described in
%the paper by Syed Alam Abbas
%====================================================================
% This function performs a Sherical ADJOINT transform on a 3D signal X
% given on the sherical coordinate system. If X is K x M x (N+1), the output will have
% (N+1) x (N+1) x (N+1) values.
% The algorithm applied here uses the Fast Fractional Fourier Transform.
%
% Synopsis: Y= Adjoint3DSphericalDFT(X)
%
% Inputs -  X      K x M x (N+1) matrix in Spherical grid,
%                  (N is assumed to be even)
%                  (M is assumed to be even)
%                  (K is assumed to be even)
% Outputs - Y      (N+1) x (N+1) x (N+1) matrix (full Cartesian grid 3D)
%
% Example:
%
%   The following is a way to verify that this works as an adjoint -
%   choosing random X and Y, one must obtain that
% The adjoint of a transformation is a unique transformation T* so that
% <Tx,y> = <x, T*y> for every x and y  where T* is an adjoint of operator T
% <y,Ax> = <x,adj(A)y>'
%
%    N=16;                                            % Always even
%    M = N +2 ;                                       % Always even
%    K = N+2;
%    X=zeros(N+1, N+1,N+1)+sqrt(-1)*zeros(N+1,N+1, N+1);   % Complex matrix of Polar Grid
%    X(1,1,:)= rand(N+1,1)+1i* rand(N+1,1);
%    %X=randn(N+1, N+1,N+1)+sqrt(-1)*randn(N+1,N+1, N+1);   % Complex matrix of Polar Grid
%    Y=zeros(K, M, N+1)+sqrt(-1)*zeros(K, M, N+1);           % Complex matrix
%    Y(1,1,:)= rand(N+1,1)+1i* rand(N+1,1);
%    %Y=randn(K, M, N+1)+sqrt(-1)*randn(K, M, N+1);           % Complex matrix
%    Ax=VectorizedCompute3DSphericalDFT(  X, K, M, zeros(2,2) );      % output K x M x N+1 Spherical polar grid
%    AHy=Adjoint3DSphericalDFT(Y);                     % Output  N+1 x N+1 x N+1   adjoint image
%    prod1 = Y.*Ax;
%    prod2 = X.*AHy;                             % ??
%    disp(abs( sum(sum(sum(prod1))) -conj(sum(sum(sum(prod2))))));
% Written on October 7th, 2015 by Syed Alam Abbas.
% Updated on October 30th, 2015 by Syed Alam Abbas
% Fixed the special Line Problem on December 7th, 2015 by Syed Alam Abbas
% Fixed a bug on December 9th, 2015 regarding switch of X-Y indexes in the exponent
%====================================================================

%% Few Initializations of required variables
global N;
global d_NoOfElements;
global lineSpacing;

[ K, M, d_NoOfElements ] = size(SphericalGrid);
N = d_NoOfElements -1;
lineSpacing = -N / 2: N / 2;
ZeroValues = zeros (d_NoOfElements, d_NoOfElements, d_NoOfElements); % Supposed to be Complex Adjoint image
Image3DAdjoint = ZeroValues ;
% FrFT1D_Uniform_Image3D;                        % 3D Image obtained after passing 1st stage of FrFT, all columns have now been operated on
% Slice2D, Slice2D_Conj;                         % 2D Slices obtained at the 2nd stage of FrFT
% Line1, Line1_Conj, Line2, Line2_Conj;		     % 1D Lines at the end of the 3rd stage of FrFT, which are final
% alpha_factor, beta_factor, gamma_factor;       % Scaling factors in X-axis, Y-axis and Z-axis respectively, which change depending on the computation of the block

%% Computing number of levels
% No of levels theta                           % These computed on fly
P = (K-2)/4;                                   % Number of levels
if(rem(K-2,4) ~= 0)
    P = ceil (P);                              % use + 1 to compute for 45 but skip for other dimension
end

% No of levels Phi
Q = (M-2)/4;                                   % Number of levels
if(rem(M-2,4) ~= 0)
    Q = ceil (Q);                              % use + 1 to compute for 45 but skip for other dimension
end

debug = 0;              % set this flag to debug

%% Processing the Reversing of the  forward Spherical Fourier Transform operation step by step
for q=1:Q  % For each level of Polar slices
    anglePhi = q*180/M;
    %% X oriented pair of Polar slices , (1+q, K+1-q)
    for p = 1:P
        angleTheta = p*180/K;
        %% XX block  -- Concentric rectangles in YZ tiled along X -axis  (1+p, M+1-p)
        alpha_factor = cosd(angleTheta) * cosd(anglePhi);											         % Scaling needed as defined in the paper
        beta_factor  = cosd(angleTheta) * sind(anglePhi);
        gamma_factor = sind(angleTheta);
        
        Line1      = SphericalGrid(p+1, q+1, :);                                %% Reversing all operations as we did in forward form
        Line2_Conj = conj(SphericalGrid(p+1, M+1 - q, :));
        Line1_Conj = conj(SphericalGrid(K+1 - p, q+1, :));						%% 1 More Swap to match !
        Line2      = SphericalGrid(K+1 - p, M+1 - q, :);
        
        %         %% Degub lines should remain commented
        %         Line1 = Line1_Conj;
        %         Line1 = zeros(size(Line1));
        %         Line2_Conj = zeros(size(Line1));
        %         %         Line1_Conj = zeros(size(Line1));
        %         Line2      = zeros(size(Line1));
        
        %% Continued processing
        % Z axis differential scaling
        Slice2D      = ComputeAdjoint_2ComplementaryLines(Line1, Line1_Conj, gamma_factor);
        Slice2D_Conj = ComputeAdjoint_2ComplementaryLines(Line2, Line2_Conj, gamma_factor);
        
        % Y axis differential scaling
        FrFT1D_Image_3D = ComputeAdjoint_2ComplementaryPlanes(Slice2D, Slice2D_Conj, beta_factor);
        
        FrFT1D_Image_3D = permute(FrFT1D_Image_3D, [2 1 3]);    %% To permute or not to permute ??? PERMUTE
        
        if (q == 1 && p == 1)                                     %% Computing Polar Slice at anglephi = 90, only needs to be computed once
            Polar2D = SphericalGrid(:, M / 2+1, :);
            properOriented2DSlice =  Adjoint2DPolarDFT( squeeze(Polar2D) );       %% Using 2D adjoint operation
            CubeSlice2D      = zeros(d_NoOfElements,d_NoOfElements,d_NoOfElements);
            for h = 1:d_NoOfElements
                CubeSlice2D(h,:,:)= properOriented2DSlice.';
            end
            CubeSlice2D_First  = permute(CubeSlice2D, [2 3 1]);
            FrFT1D_Uniform_Image3D = permute(CubeSlice2D_First, [2 3 1]);
            FrFT1D_Uniform_Image3D = permute (FrFT1D_Uniform_Image3D, [ 2 1 3]);
            Image3DAdjoint = Image3DAdjoint + FrFT1D_Uniform_Image3D;                          % collecting contribution from every level
        end
        
        % X axis final uniform scaling
        FrFT1D_Uniform_Image3D = VectorizedFrFT_Centered(FrFT1D_Image_3D, -alpha_factor); % Observe the sign change from -ve to +ve
        %                 FrFT1D_Uniform_Image3D = permute (FrFT1D_Uniform_Image3D , [2 1 3]);
        Image3DAdjoint = Image3DAdjoint + FrFT1D_Uniform_Image3D;                          % collecting contribution from every level
        
        %% If theta = 45, skip the next block which is redundant
        if (angleTheta == 45 )
            continue;
        end
        %% XZ block  -- Concentric rectangles in
        alpha_factor = sind(angleTheta) * cosd(anglePhi);		 									         % Scaling needed as defined in the paper
        beta_factor  = sind(angleTheta) * sind(anglePhi);
        gamma_factor = cosd(angleTheta);
        
        Line1      = SphericalGrid(K/2+1 - p, 1+q, :);
        Line1_Conj = SphericalGrid(K/2+1 - p, M+1 - q, :);
        Line2_Conj = SphericalGrid(K/2+1 + p, q+1, :);
        Line2      = SphericalGrid(K/2+1 + p, M+1 - q, :);
        
        %         Line1 = Line2;          %% Degub lines should remain commented
        %         Line1 = zeros(size(Line1));
        %         Line2_Conj = zeros(size(Line1));
        %         Line1_Conj = zeros(size(Line1));
        %         Line2 = zeros(size(Line1));
        
        
        % X axis differential scaling
        Slice2D = ComputeAdjoint_2ComplementaryLines(Line1, flipud(conj( squeeze(Line1_Conj))), alpha_factor);
        Slice2D_Conj = ComputeAdjoint_2ComplementaryLines(Line2,flipud(conj( squeeze(Line2_Conj))),alpha_factor);
        
        % Y axis differential scaling
        FrFT1D_Image_3D = ComputeAdjoint_2ComplementaryPlanes(Slice2D, Slice2D_Conj, beta_factor);
        FrFT1D_Image_3D  = permute (FrFT1D_Image_3D , [2 1 3]);
        
        if (q == 1 && p == 1)         %% Computing One special Polar Slice at angletheta = 0 , only need to be computed once
            Polar2D = SphericalGrid(1, :, :);
            Polar2D(1,M/2+1,:)= zeros(1,d_NoOfElements);              % Not repeating pattern line of SphericalGrid(1, M/2+1, :) since it was used in the previous XX block
            properOriented2DSlice =  Adjoint2DPolarDFT(squeeze(Polar2D));
            CubeSlice2D      = zeros(d_NoOfElements,d_NoOfElements,d_NoOfElements);
            for h = 1:d_NoOfElements
                CubeSlice2D(h,:,:)= properOriented2DSlice.';
            end
            CubeSlice2D_First  = permute(CubeSlice2D, [3 1 2]);
            FrFT1D_Uniform_Image3D = permute(CubeSlice2D_First, [3 1 2]);
            FrFT1D_Uniform_Image3D = permute (FrFT1D_Uniform_Image3D, [2 1 3]);               % Just make it work , then figure out how to fix it simply what is the permutation
            Image3DAdjoint = Image3DAdjoint + FrFT1D_Uniform_Image3D;                          % collecting contribution from every level
        end
        
        
        
        %% If Phi = 45, skip the next block which is redundant
        if (anglePhi == 45)   % If angle is 45 we have already computed it above
            % Z axis final uniform scaling alone for the upper grid XZ block else
            % compute it later after the YZ block processing since they
            % share gamma uniform scaling in the Z-axis
            % Z axis final uniform scaling
            FrFT1D_Uniform_Image3D = VectorizedFrFT_Centered(FrFT1D_Image_3D, -gamma_factor); % Observe the sign change from -ve to +ve
            FrFT1D_Uniform_Image3D = permute (FrFT1D_Uniform_Image3D , [3 2 1]);
            %         if(p == 1 && q == 1 )
            %         if(p == 2 && q == 1)
            %           if(p == 2 && q == 2)
            %           if(p == 1 && q == 2)
            %         if((p == 2 || p == 1) && (q == 1 || q == 2))
            
            Image3DAdjoint = Image3DAdjoint + FrFT1D_Uniform_Image3D;                          % collecting contribution from every level
            %         end
            continue;
        end
        %% YZ block -- Concentric rectangles in
        alpha_factor = sind(angleTheta) * sind(anglePhi);		   % Scaling needed as defined in the paper
        beta_factor  = sind(angleTheta) * cosd(anglePhi);
        gamma_factor = cosd(angleTheta);                             % see the same scale so use the previously computed grid FrFT 1D
        
        Line1      = SphericalGrid(K / 2+1 - p, M / 2+1 - q, :) ;
        Line1_Conj = SphericalGrid(K / 2+1 - p, M / 2+1 + q, :);
        Line2_Conj = SphericalGrid(K/ 2+1 + p, M/ 2+1 - q, :);
        Line2      = SphericalGrid(K/ 2+1 + p, M / 2+1 + q, :);
        
        % X axis differential scaling
        Slice2D = ComputeAdjoint_2ComplementaryLines(Line1, flipud(conj( squeeze(Line1_Conj))), alpha_factor);
        Slice2D_Conj = ComputeAdjoint_2ComplementaryLines(Line2, flipud(conj( squeeze(Line2_Conj))),alpha_factor);
        
        % Y axis differential scaling and addition from XZ block
        FrFT1D_Image_3D_YZ =  ComputeAdjoint_2ComplementaryPlanes(Slice2D, Slice2D_Conj, beta_factor);
        FrFT1D_Image_3D_YZ = permute (FrFT1D_Image_3D_YZ , [2 1 3]);
        FrFT1D_Image_3D = FrFT1D_Image_3D + FrFT1D_Image_3D_YZ ;                                     %% Adding contribution form XZ block
        
        % Z axis final uniform scaling - combination of YZ and XZ blocks
        FrFT1D_Uniform_Image3D = VectorizedFrFT_Centered(FrFT1D_Image_3D, -gamma_factor); % Observe the sign change from -ve to +ve
        FrFT1D_Uniform_Image3D = permute (FrFT1D_Uniform_Image3D , [3 2 1]);
        Image3DAdjoint = Image3DAdjoint + FrFT1D_Uniform_Image3D;                          % collecting contribution from every level
        
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
        beta_factor  = cosd(angleTheta) * cosd(anglePhi);
        gamma_factor = sind(angleTheta);
        
        Line1      = SphericalGrid(p+1, M/ 2+1 - q, :);
        Line2      = SphericalGrid(p+1, M/ 2+1 + q, :) ;
        Line1_Conj = conj(SphericalGrid(K+1 - p, M/ 2+1 - q, :));
        Line2_Conj = conj(SphericalGrid(K+1 - p, M/ 2+1 + q, :));
        
        % Z axis differential scaling
        Slice2D = ComputeAdjoint_2ComplementaryLines(Line1, Line1_Conj, gamma_factor);
        Slice2D_Conj = ComputeAdjoint_2ComplementaryLines(Line2, Line2_Conj,gamma_factor);
        
        % X axis differential scaling
        FrFT1D_Image_3D = ComputeAdjoint_2ComplementaryPlanes(Slice2D, Slice2D_Conj, alpha_factor);
        
        if(p == 1 && q == 1)% Computing One special Polar Slice at  anglephi = 0 , only need to be computed once
            Polar2D = SphericalGrid(:, 1, :) ;
            Polar2D(1,1,:)= zeros(1,d_NoOfElements);              % Not repeating pattern line of SphericalGrid(1, 1, :) since it was used in the previous XZ block
            properOriented2DSlice =  Adjoint2DPolarDFT( squeeze(Polar2D) );       %% Using 2D adjoint operation
            CubeSlice2D      = zeros(d_NoOfElements,d_NoOfElements,d_NoOfElements);
            for h = 1:d_NoOfElements
                CubeSlice2D(h,:,:)= properOriented2DSlice;
            end
            CubeSlice2D = permute (CubeSlice2D, [ 2 1 3]);
            Image3DAdjoint = Image3DAdjoint + CubeSlice2D;                          % collecting contribution from every level
            
            %% Special operation was missing  earlier now in place
            Polar2D = SphericalGrid(K/2+1, :, :) ;
            Polar2D = squeeze(Polar2D);
            Polar2D(1,:)= zeros(1,d_NoOfElements);              % Not repeating pattern line of SphericalGrid(1, 1, :)
            Polar2D(M/2+1,:)= zeros(1,d_NoOfElements);           % Not repeating pattern line of SphericalGrid(1, M/2+1, :) since it was used in the previous XX block
            
            %             for d = 1:N+1            % Should we optimize this , YES (see below), M x (N+1) x (N+1) complexity
            %                 for m = 1:M
            %                     for n=1:N+1
            %                         SpecialLineZ(d) = SpecialLineZ(d) + Polar2D(m,n)*exp(+1i*2*pi*(lineSpacing(d)*lineSpacing(n))/(N+1));
            %                     end
            %                 end
            %             end
            
            SpecialLineZ= zeros(1,N+1);
            for d = 1:N+1            % Should we optimize this ???? M x (N+1) x (N+1) complexity
                SpecialLineZ(d) = SpecialLineZ(d) + sum(Polar2D*exp(+1i*2*pi*(lineSpacing(d)*lineSpacing')/(N+1)));
            end
            
            properOriented2DSlice = repmat(squeeze(SpecialLineZ), N+1,1);
            CubeSlice2D      = zeros(d_NoOfElements,d_NoOfElements,d_NoOfElements);
            for h = 1:d_NoOfElements
                CubeSlice2D(h,:,:)= properOriented2DSlice;
            end
            Image3DAdjoint = Image3DAdjoint + CubeSlice2D;                          % collecting contribution from every level
        end
        % Y-axis final uniform scaling
        FrFT1D_Image_3D = permute (FrFT1D_Image_3D, [ 2 1 3]);
        FrFT1D_Uniform_Image3D = VectorizedFrFT_Centered(FrFT1D_Image_3D, -beta_factor);      % Observe the sign change from -ve to +ve
        FrFT1D_Uniform_Image3D = permute (FrFT1D_Uniform_Image3D, [ 2 1 3]);
        Image3DAdjoint = Image3DAdjoint + FrFT1D_Uniform_Image3D;                          % collecting contribution from every level
    end
    
end

% Image3DAdjoint = DebugAdjoint_3DDFT;
if(debug)
    disp('Full debug information')
    squeeze(Image3DAdjoint(1,:,:))           % Along First Axis
    squeeze(DebugAdjoint_3DDFT(1,:,:))
    
    squeeze(Image3DAdjoint(:,2,:))           % Along Second Axis
    squeeze(DebugAdjoint_3DDFT(:,2,:))
    
    squeeze(Image3DAdjoint(:,:,3))           % Along Third Axis
    squeeze(DebugAdjoint_3DDFT(:,:,3))
    
    squeeze(Image3DAdjoint - DebugAdjoint_3DDFT)
end

    function [FrFT1D_Image_3D] = ComputeAdjoint_2ComplementaryPlanes(Slice2D, Slice2D_Conj, ColumnScale)
        
        lineSpacing_Cube_TiledLevel = repmat(lineSpacing'*lineSpacing,[1 1 d_NoOfElements]);
        Column_ScaleFactor = exp (+2*1i * pi * lineSpacing_Cube_TiledLevel * ColumnScale / d_NoOfElements);   % Observe the sign change from -ve to +ve
        
        CubeSlice2D_first = repmat(Slice2D, [1 1 d_NoOfElements]);
        CubeSlice2DConj  = repmat(Slice2D_Conj, [1 1 d_NoOfElements]);
        
        CubeSlice2D_first =  permute (CubeSlice2D_first, [ 3 2 1]);
        CubeSlice2DConj   =  permute (CubeSlice2DConj, [ 3 2 1]);
        
        FrFT3DBlock1 = Column_ScaleFactor .* CubeSlice2D_first;
        FrFT3DBlock2 = conj(Column_ScaleFactor) .* CubeSlice2DConj;
        
        FrFT1D_Image_3D = FrFT3DBlock1 + FrFT3DBlock2;
    end

    function [FrFTFrFT_Image2D] = ComputeAdjoint_2ComplementaryLines(Line1D, Line1D_Conj, ColumnScale)
        lineSpacing_Square_TiledLevel = lineSpacing'*lineSpacing;
        Column_ScaleFactor = exp(+2*1i * pi * lineSpacing_Square_TiledLevel * ColumnScale / d_NoOfElements);      % Observe the sign change from -ve to +ve
        
        tiledLine1D      = squeeze(Line1D) * ones(1, d_NoOfElements) ;
        tiledLine1D_Conj = flipud(conj(squeeze(Line1D_Conj))) * ones(1, d_NoOfElements) ;        % This thing is needed, check the solution
        
        FrFT_2DBlock1 = Column_ScaleFactor .* tiledLine1D.';
        FrFT_2DBlock2 = conj(Column_ScaleFactor) .* tiledLine1D_Conj.';
        
        FrFTFrFT_Image2D = FrFT_2DBlock1 + FrFT_2DBlock2;
    end
end

