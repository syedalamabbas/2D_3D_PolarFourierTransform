function [ Image_Final ] = Specialized2DInverse( DFT_Polar,DFT_Corners, N, M )
%====================================================================
% SPECIALIZED2DINVERSE Takes in the input DFT information from polar grid
% and the cartesian grid with just the corners
% DFT_Polar   = M x (N+1) grid of complex numbers
% DFT_Corners = (N+1) x (N+1) with only central elements as zero
% N+1         - is the image size desired
% M           - is the number of angles
% Written on June 29th, 2015 by Syed Alam Abbas.
% Update on December 13th, 2015 by Syed Alam Abbas, functionality on BTTB to BCCB functions 
%====================================================================
% This function actually needs improvement and better functionality
% implementation wise, maybe follow the inverse techinique used in NFFT
% paper 2007 Applied Computational Harmonic Analysis
if(nargin == 1),
    DFT_Corners = zeros(N+1);     % Assuming nulled corners
end

persistent BCCB_EigenValues;

%% Polar Adjoint Comutation for System 1
ImpulseImage = zeros(N+1);
ImpulseImage(1,1) = 1;

W         = sqrt(abs(-N/2:N/2)/2)/(N+1);    % Fourier based preconditioner from Amir's 2001 paper Fast Slant Slack
W(N/2+1)  = sqrt(1/8)/(N+1);
W         = repmat(W, M,1);

figure, imagesc(W)
title('Weighting matrix')
ylabel(['Angles, M= ',num2str(M)] )
xlabel(['Radial line, N+1 = ',num2str(N+1)])
colorbar 
% Const_factor=M*(((N+1)/2.0)*((N+1)/2.0)+1.0/4.0);
% W = abs(-N/2:N/2)/Const_factor;
% W(N/2+1) = 1.0/4.0/Const_factor;
% W = ones(M,1)*W;
% W = ones(M,N+1);     % No weights check??

% oneDWeights = abs(linspace(-1/2,1/2,(N+1))); % This is a little involved dont use it //
% W         = diag(pi*oneDWeights)/(M*(N+1));        % Weights defined by 2007 paper , this is j
% PolarPartImpulse =  Adjoint2DPolarDFT(W .* Compute2DPolarDFT( ImpulseImage,M) );

PolarPartImpulse =  Adjoint2DPolarDFT(W.*W.* Compute2DPolarDFT( ImpulseImage, M ) );  % Use this for Amir's weights

%% Corner Cartesian with Polar Nulling Comutation for System 2
I_DFT_Corner = fftshift(fft2(ImpulseImage));
GridSpacing = -N/2:N/2;
[CartesianGridX, CartesianGridY]= meshgrid(GridSpacing,GridSpacing);
Weight=(CartesianGridX.^2 + CartesianGridY.^2 < (N/2)^2);                  % Polar centered points nulling
Pos = find(Weight(:));
zeroArray = zeros(length(Pos),1);
I_DFT_Corner(Pos) = zeroArray;
figure, imagesc(abs(I_DFT_Corner))
CornerPartImpulse = ifft2(ifftshift(I_DFT_Corner));      % This is not necessarily all zeros, so must keep it
figure, imagesc(CornerPartImpulse)
%% Computing the collective contribution ---- Polar (System 1) + Corner (System 2) and making the larger BCCB matrix
BTTB_ColumnsIn2DShape = PolarPartImpulse + CornerPartImpulse;
BCCB_ColumnIn2DShape = makeBCCBFromBTTB_2Dextension(BTTB_ColumnsIn2DShape,N);
BCCB_EigenValues = fft2(BCCB_ColumnIn2DShape);       % For a BCCB matrix  we have the spectral decomposition, A = F* Eig F, F* is the ifft2, F is fft2

%% Forming the actual solution  A x = b
% PolarPart = Adjoint2DPolarDFT(W .*  DFT_Polar);      % A^H_p using Markus's weights

PolarPart = Adjoint2DPolarDFT(W.*W.* DFT_Polar);      % A^H_p  using Amir's weights again
CornersPart =  ifft2(ifftshift(DFT_Corners));         % A^H_cc 
CornersPart = zeros(size(PolarPart));                 % Setting default value for corners

sum2D_b = PolarPart + CornersPart ;                   % b_2D = A^H_p + A^H_cc    
b = reshape(sum2D_b, [(N+1)^2,1]);

stackedImageColumnwise = cgs(@A,b,10^-4,2000);               % MATLAB function solution
% stackedImageColumnwise = cgsolve(@A,b,10^-4,200,1);             % Verbose solution using custom program using different iterations
Image_Final = reshape ( stackedImageColumnwise, [N+1, N+1] );
Image_Final = real(Image_Final);
Image_Final = fliplr(rot90(Image_Final,-1));                 % Dont know why this rotation is required, but it works ! 
    
    
    function b_iterate_1D = A( x_iterate_1D )                     % A = (A^H_p A_p + A_cc A^H_cc)
      b_iterate_1D =  BTTB_2D_Matrix_Vector(x_iterate_1D,BCCB_EigenValues,N );
    end
end

% function b_iterate_1D = BTTB_2D_Matrix_Vector(x_iterate_1D,BCCB_EigenValues,N )
%     inputImage       = reshape ( x_iterate_1D, [N+1, N+1]);           % Column stacking unfolding
%     inputImage       = padarray(inputImage, [N-1,N-1], 'post');
%     b_iterate2D      = ifft2( BCCB_EigenValues .* fft2(inputImage) ); % Matrix computations with the BCCB matrix multiplication
%     b_iterate2D      = real( b_iterate2D );                           % Due to round off error the computed solutions may contain complex parts
%     b_iterateRefined = b_iterate2D( 1:N+1,1:N+1 );                    % Taking just the first part of the BCCB matrix multiplication leaving last rows and cols
%     b_iterate_1D     = reshape( b_iterateRefined, [(N+1)^2, 1] );     % Column stacking
% end

% function BCCB_ColumnIn2DShape = makeBCCBFromBTTB_2Dextension(BTTB_ColumnsIn2DShape,N)           % N is even and same definition
% % This function assumes that the BTTB is symmetric real matrix
% BCCB_ColumnIn2DShape = zeros(2*N, 2*N);                                 % Making a full block circulant with circulant blocks
% for k= 1:N+1
%     bttb_col =  BTTB_ColumnsIn2DShape(:,k);
%     flippedColumn = flipud(bttb_col);
%     BCCB_ColumnIn2DShape(:,k) = [bttb_col; flippedColumn(2:N)];         % periodic extension
% end
% flippedIndexes = fliplr(2:N);                    % Making full matrix block circulant
% for k= 1:N-1
%     bttb_col = BTTB_ColumnsIn2DShape(:,flippedIndexes(k));
%     flippedColumn = flipud(bttb_col);
%     BCCB_ColumnIn2DShape(:,k+N+1) = [bttb_col; flippedColumn(2:N)];     % periodic extension
% end
% end


% %% Testing the two structure like A = A_1 : A_2, is the augumented matrix which has two systems 1 and systems 2, 
% A^H * A, Combined system
% A^H_1 * A_1 + A^H_2 * A_2, Summing the two systems  and checking if it matches the solution  
% N = 8;
% a= sym('a',[1,N]);
% b= sym('b',[1,N]); 
% c= sym('c',[1,N]); 
% d= sym('d',[1,N]);
% e= sym('e',[1,N]); 
% f= sym('f',[1,N]);
% g= sym('g',[1,N]);
% h= sym('h',[1,N]);
% i = sym('i',[1,N]);
% j = sym('j',[1,N]);
% 
% A_1 = [a;b;c;d;e]
% A_2 = [f;g;h;i;j]
% A = [a;b;c;d;e;f;g;h;i;j]
% A'*A - (A_1'*A_1 + A_2'*A_2)   % The error to be zero