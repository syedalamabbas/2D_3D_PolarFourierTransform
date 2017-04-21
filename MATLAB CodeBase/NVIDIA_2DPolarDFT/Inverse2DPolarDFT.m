function [ ImageFinal ] = Inverse2DPolarDFT( Polar_Grid, halfSizeDesired, accuracy )

% % PolarGridCorners1,PolarGridCorners2, PolarGridCorners3,PolarGridCorners4,
%INVERSE2DPOLARDFT is the inverse polar grid computation, given Polar_grid
%of size M x (N+1), where M is the number of angles (even)
% N+1 is the number of elements in each anglular slice of the 2D polar grid
% Solution is based on Fourier domain interpolation on circles followed by
% FrFT in radial and angular directions

if nargin==2,
    accuracy=1e-5;   % Default accuracy
end;

[M, sizeN_plus1 ] = size(Polar_Grid);
N = sizeN_plus1 -1;
ImageFinal = zeros(N+1,N+1);

% %% T Chan's preconditioner Calculation, only dependent on A^H A
% ImpulseImage = zeros(N+1,N+1);
% ImpulseImage(1,1) = 1;            % Making it an impulse or a point spread function, 
%   % Special unique column of the BTTB matrix which is formed by product (A^H * A). This can be computed only once at the beginning for different inputs
%  BTTB_ImpluseColumnsIn2DShape = Adjoint2DPolarDFT( Compute2DPolarDFT( ImpulseImage,  M ));     % It is (N+1) x (N+1) size, it is a impluse response of the linear system
%  BCCB_ImpulseColumnIn2DShape = ComputeBCCBPreConditioner( BTTB_ImpluseColumnsIn2DShape );
  
% %% Highly Specialized Direct Solution
% ImageFinal(1,1) = 1;            % Making it an impulse or a point spread function, 
% % Special unique column of the BTTB matrix which is formed by product (A^H * A). This can be computed only once at the beginning for different inputs
% BTTB_ColumnsIn2DShape = Adjoint2DPolarDFT( Compute2DPolarDFT( ImageFinal,  M ));     % It is (N+1) x (N+1) size, it is a impluse response of the linear system
% figure, imagesc(real(BTTB_ColumnsIn2DShape))
% title('Impulse response or point Spread function')
% % Transforming the first column of BTTB into BCCB column, requires 
% % reflexive extentions per block as well as for entire matrix
% BCCB_ColumnIn2DShape = makeBCCBFromBTTB_extension(BTTB_ColumnsIn2DShape,N);
% figure, imagesc(real(BCCB_ColumnIn2DShape))
% title('Block Circulant (BCCB) extension')
% 
% %% Solving x = (A^H*A)^-1 A^H y   the full linear system like x = A^-1 b
% b = Adjoint2DPolarDFT( Polar_Grid);   % (N+1) x (N+1)   A^H y
% padded_b =  padarray(b, [N-1 N-1], 'post');                          % extension of the input
% Eigen_values_BCCB = fft2(BCCB_ColumnIn2DShape);      
% ImageFullFinal = ifft2( fft2(padded_b) ./  Eigen_values_BCCB);
% figure, imagesc( ( real(ImageFullFinal)))
% ImageFinal  = real(ImageFullFinal(1:N+1, 1:N+1));   % We use this real() since due to overflow error small imaginary numbers may surface




%% preconditioner 1
fullSize  = sizeN_plus1;
W         = sqrt(abs(-N/2:N/2)/2)/fullSize;    % Fourier based preconditioner from Amir's paper
% % W         = sqrt(2.5*abs(-N/2:N/2))/fullSize; 
W(N/2+1)  =  sqrt(1/8)/fullSize;
% W(N/2+1)  = sqrt(1/1)/fullSize;
 
% W = pi/M* [-N/2:N/2]/ (N+1);      % New Analytical weights from Markus's paper
% W(N/2+1) = pi/4*(1/(M*(N+1)));
 W        = ones(M,1)*W ;
  
 
 
 
 %% preconditioner 2 most generic , inverse of the diagonals
%  W = 1/ (M*(N+1))*eye(N+1);    % Generic Preconditioner ???

%% Simple Gradient Descent
Delta=1;
count=0;
maxIterations = 10;
while Delta>accuracy && count < maxIterations, 
    
    %% No preconditioner 
    
%      Err =  (Compute2DPolarDFT( ImageFinal,  M ) - Polar_Grid);
%   D = Adjoint2DPolarDFT(Err).';
%%   Using preconditioner 1
  Err = W.* (Compute2DPolarDFT( ImageFinal,  M ) - Polar_Grid);
  
%   UniformSpacing1 = [-halfSizeDesired halfSizeDesired];
%   deltaTheta = 180/M;     % Angular sampling rate
%   angles = [0:deltaTheta:180-deltaTheta];
%   for k = 1:length(angles)
%       angle = angles(k);
%       if (  angle <= 45 || angle > 135 )
%           newGridSpacing = UniformSpacing1 ./ cosd(angle);       % BH
%           if(angle > 135)
%               halfSpacing = 1:1:newGridSpacing(1);
%           else
%               halfSpacing = 1:1:newGridSpacing(end);
%           end
%       else if ( angle > 45 || angle <= 135 )
%               newGridSpacing = UniformSpacing1 ./ sind(angle);    %BV
%               halfSpacing = 1:1:newGridSpacing(end);
%           end
%       end
%       
%       
%       %  newGridSpacing = [-fliplr(halfSpacing) 0 halfSpacing]; <=   UniformSpacing2
%       
%       CurrentRadialLine = Err(k,:);                                                  % Retrieving
%       indexesToNull = [ (-[halfSpacing(end)+1:N/2]+ N/2+1)  ([halfSpacing(end)+1:N/2]+ N/2+1)  ] ;       % Processing
%       CurrentRadialLine(indexesToNull) = zeros(1,length(indexesToNull));                    % Updating
%       Err(k,:) =  CurrentRadialLine ;                                                 % Replacing
%   end
  
  D = Adjoint2DPolarDFT(W.* Err).';

  %% Discarding all information beyond the image seen and desired central elements only
%    startIndex = N/2 - halfSizeDesired +1;
%    endIndex = startIndex + (2*halfSizeDesired);
%    centralIndexes = [1:startIndex-1,endIndex+1:fullSize];
%    D(centralIndexes,:)= zeros(length(centralIndexes),fullSize);   % Removing all info beyond central elements
%    D(:, centralIndexes) = zeros(fullSize,length(centralIndexes));
%   Mean = 0;          % Very Suceptible to noise due to 
%   Var = .0001;
%   RealPart = imnoise(real(Polar_Grid),'gaussian',Mean,Var);
%   ImagPart = imnoise(imag(Polar_Grid),'gaussian',Mean,Var);
%   Polar_Grid = RealPart + 1i* ImagPart; 
  
% retImage = Compute2DPolarDFT( ImageFinal,  M );
%   Err =  ( retImage- Polar_Grid); 
%   
%   [CurrentPolarCorners,WeightsCorners,~] = Compute2DPolarCornersDFT( ImageFinal,  M );
%   Err1 = CurrentPolarCorners - PolarGridCorners1;
%   mat1 = Adjoint2DPolarCornersDFT( WeightsCorners' .*WeightsCorners' .* Err1, N,M ).';
%   
%   [CurrentPolarCorners,WeightsCorners,~] = Compute2DPolarCornersDFT( flipud(ImageFinal),  M );
%   Err2 = CurrentPolarCorners - PolarGridCorners2;
%   mat2 = Adjoint2DPolarCornersDFT( WeightsCorners' .*WeightsCorners' .* Err2, N,M ).';
%   
%   [CurrentPolarCorners,WeightsCorners,~] = Compute2DPolarCornersDFT( fliplr(ImageFinal),  M );
%   Err3 = CurrentPolarCorners - PolarGridCorners3;
%   mat3 = Adjoint2DPolarCornersDFT( WeightsCorners' .*WeightsCorners' .* Err3, N,M ).';
%   
%   [CurrentPolarCorners,WeightsCorners,~] = Compute2DPolarCornersDFT(  rot90(ImageFinal,2),  M );
%   Err4 = CurrentPolarCorners - PolarGridCorners4;
%   mat4 = Adjoint2DPolarCornersDFT( WeightsCorners' .*WeightsCorners' .* Err4, N,M ).';
%   
% %   D = ( Adjoint2DPolarDFT(W.*W.*  Err)).';
%   
%   D = ( Adjoint2DPolarDFT(W.*W.*  Err)).' + mat1 + mat2+ mat3+ mat4;
%   D = ( Adjoint2DPolarDFT( Err)).' + mat1;
  
  

  
%   D = ( W.*W.* Adjoint2DPolarDFT( Compute2DPolarDFT( ImageFinal,  M ) - Polar_Grid)).';
%   
    
%   %% Using T. Chan's BCCB preconditioner
%   S = fft2(BCCB_ImpulseColumnIn2DShape);  % Filter
%   
%   X1 = Adjoint2DPolarDFT( Polar_Grid);            % P* A^H (b)
%   D1 = ifft2 ( S .* fft2 (X1));                        % Multiplication by a Preconditioner BCCB matrix, Fast solution ----> BCCB_Matrix * A^H (Ax - b)
%   D1 = real(D1);
%   
%   X = Adjoint2DPolarDFT( Compute2DPolarDFT( ImageFinal,  M ));            % P* A^H (Ax)
%   D2 = ifft2 ( S .* fft2 (X));                        % Multiplication by a Preconditioner BCCB matrix, Fast solution ----> BCCB_Matrix * A^H (Ax - b)
%   D2 = real(D2);
%   
%   D = D2 - D1;
  
%%   Using preconditioner 2
%     D = W* (Adjoint2DPolarDFT((Compute2DPolarDFT( ImageFinal,  M ) - Polar_Grid)).');

%% Computing the solution
%     Delta=norm(D(startIndex:endIndex,startIndex:endIndex));
    Delta=norm(D);
    disp(['At iteration ',num2str(count),' the difference matrix norm is ',num2str(Delta)]);
    mu=1/sizeN_plus1;
    Temp=ImageFinal-mu*D; 
    ImageFinal=Temp;
    count=count+1;  
    if(count >= 1 && count <= 3)
        figure, imshow(real(ImageFinal), [])
        str = strcat('Retrieval at iteration no.' , num2str(count));
        title (str );
      
    end
end; 
disp(['Number of required iterations is ',num2str(count)]);

figure, imshow(real(ImageFinal), [])
% figure, imagesc(imadjust(real(ImageFinal)))
str = strcat('Retrieval at iteration no.' , num2str(count));
title (str );
 
 
%% GMRES Solution  
% StackedRadialLines = reshape ( W* Polar_Grid (1:sizeN_plus1, :), [sizeN_plus1^2,1]);
% b = (  StackedRadialLines );
% 
% stackedImageColumnwise = gmres(@afun,b,sizeN_plus1^2,10^-13,sizeN_plus1^2);          % using gmres
% 
% % stackedImageColumnwise = cgs(@afun,b,10^-13,sizeN_plus1);          % using MATLAB cg
% 
% ImageFinal = reshape ( stackedImageColumnwise, [sizeN_plus1, sizeN_plus1]);
% 
%   
%   function y = afun(x)
%           % sizeN_plus1 Must be an integer and must be odd
%          inputImage = reshape ( x, [sizeN_plus1, sizeN_plus1]);         % Column stacking
%          PolarGridIterate = VectorizedCompute2DPolarDFT( inputImage,  sizeN_plus1 +1 );
%          localStackedRadialLines  = reshape ( W * PolarGridIterate (1:sizeN_plus1, :) ,[ sizeN_plus1^2,1]);
%          y = (localStackedRadialLines);
%   end

end 

function AppendedTwoDStructure = makeBCCBFromBTTB_extension (TwoDStructure,N)
% This function assumes that the BTTB is symmetric real matrix
BCCB_Column = [];                                 % Making a full block circulant with circulant blocks
for k= 1:N+1
    bttb_col =  TwoDStructure(:,k);
    flippedColumn = flipud(bttb_col);
    BCCB_Column = [BCCB_Column; bttb_col; flippedColumn(2:N) ];            % periodic extension
end

flippedIndexes = fliplr(2:N);                   % Making full matrix block circulant
for k= 1:N-1
    bttb_col = TwoDStructure(:,flippedIndexes(k));
    flippedColumn = flipud(bttb_col);
    BCCB_Column = [BCCB_Column; bttb_col; flippedColumn(2:N) ];         % periodic extension
end

AppendedTwoDStructure = reshape(BCCB_Column, [2*N 2*N ]);
end
