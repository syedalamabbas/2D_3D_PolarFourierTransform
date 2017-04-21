function [ Image3D_Final ] = Specialized3DInverse( DFT_Spherical_Polar )
%====================================================================
% SPECIALIZED2DINVERSE Takes in the input DFT information from polar grid
% and the cartesian grid with just the corners
% DFT_Spherical_Polar   = K x M x (N+1) grid of complex numbers
% N+1         - is the image size desired
% K           - is the number of angles theta this is always even
% M           - is the number of angles phi this is always even
% Image3D_Final - Final output of the Inverse operation
% Written on December 10th, 2015 by Syed Alam Abbas.
% Updated on December 13th, 2015 , BTTB to BCCB 3-level difficult composition
%====================================================================
[K,M,SizeX ]= size(DFT_Spherical_Polar);
N = SizeX-1;

%% Spherical Polar Adjoint Comutation for System 1
ImpulseImage = zeros(N+1,N+1,N+1);
ImpulseImage(1,1,1) = 1;

W         = sqrt(abs(-N/2:N/2)/2)/(N+1);    % Fourier based preconditioner from Amir's 2001 paper Fast Slant Slack
W(N/2+1)  = sqrt(1/8)/(N+1);
W         = repmat(W,M,1);                % 3D weighting matrix size K x M x (N+1)
W         = repmat(W,[1 1 K]); 
W         = permute (W, [2 3 1]); 
W         = permute (W, [2 3 1]); 
SphericalPolarPartImpulse =  Adjoint3DSphericalDFT(W.*W.* VectorizedCompute3DSphericalDFT( ImpulseImage,K, M, [] ) , []);  % Use this for Amir's weights

%% Corner Cartesian with Spherical Polar Nulling Comutation for System 2    -TODO  Need to work on this implementation , this is not correct
Impulse_DFT_Corner = fftshift(fftn(ImpulseImage));
GridSpacing = -N/2:N/2;
[CartesianGridX, CartesianGridY, CartesianGridZ]= meshgrid(GridSpacing,GridSpacing,GridSpacing);
Weight=(CartesianGridX.^2 + CartesianGridY.^2 + CartesianGridZ.^2 < (N/2)^2);                  % Polar centered points nulling
Pos = find(Weight(:));
zeroArray = zeros(length(Pos),1);
Impulse_DFT_Corner(Pos) = zeroArray;
% CornerPartImpulse = ifft2(ifftshift(Impulse_DFT_Corner));
CornerPartImpulse = zeros(size(ImpulseImage));
%% Computing the collective contribution ---- Polar (System 1) + Corner (System 2) and making the larger BCCB matrix
BTTB_ColumnsIn3DShape = SphericalPolarPartImpulse + CornerPartImpulse;
BCCB_3D_EigenValues = computeEigenValues3D(BTTB_ColumnsIn3DShape,N);      

%% Forming the actual solution  A x = b
PolarPart = Adjoint3DSphericalDFT(W.*W.* DFT_Spherical_Polar , []);      % A^H_p  using Amir's weights again
% CornersPart =  ifft2(ifftshift(DFT_Corners));         % A^H_cc 
CornersPart = zeros(size(PolarPart));

sum2D_b = PolarPart + CornersPart ;                   % b_2D = A^H_p + A^H_cc    
b = reshape(sum2D_b, [(N+1)^3,1]);
stackedImageColumnwise = cgs(@A,b,10^-4,200);               % MATLAB function solution
% stackedImageColumnwise = cgsolve(@A,b,10^-4,200,1);             % Verbose solution using custom program using different iterations
Image3D_Final = reshape ( stackedImageColumnwise, [N+1, N+1, N+1] );
Image3D_Final = real(Image3D_Final);
    function b_iterate_1D = A( x_iterate_1D )                     % full matrix multiplication A = (N+1)^3 x (N+1)^3
         b_iterate_1D =  BTTB_3D_Matrix_Vector(x_iterate_1D,BCCB_3D_EigenValues,N );
    end

end

function BCCB_3D_EigenValues = computeEigenValues3D(BTTB_ColumnsIn3DShape,N)
    BCCB_3D_EigenValues = zeros(2*N,2*N,N+1); 
    for k = 1:N+1
        BCCB_3D_EigenValues(:,:,k) = fft2(makeBCCBFromBTTB_2Dextension(BTTB_ColumnsIn3DShape(:,:,k),N));   % For a BCCB matrix  we have the spectral decomposition, A = F* Eig F, F* is the ifft2, F is fft2
    end
end
