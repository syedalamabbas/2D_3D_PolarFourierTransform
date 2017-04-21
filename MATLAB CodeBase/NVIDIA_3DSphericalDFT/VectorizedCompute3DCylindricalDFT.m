function [ CylindricalGrid ] = VectorizedCompute3DCylindricalDFT( inputImage,  noOfAnglesPhi )
%VECTORIZEDCOMPUTE2DCYLINDRICALDFT Transforms the fucntion using the prescribed
% scheme as discussed in the paper by Syed Alam Abbas
%====================================================================
% Written on October 9th, 2015 by Syed Alam Abbas.
%====================================================================
I = inputImage;
N = size(I,1)-1;
M = noOfAnglesPhi;

CylindricalGrid = zeros(M, N+1, N+1);                % Z -coordinates vs  Polar slices: No of angles vs. Radial data

%% Simple operation for computing Cylindrical grid
ReorderedImage_OperateColumns = permute(I, [3 2 1]);                                %  achieving with Swap Z-axis as columns now
FrFT1D_Uniform_Image3D = VectorizedFrFT_Centered(ReorderedImage_OperateColumns,1);  %  compute Z-axis centered FFT
FrFT1D_Uniform_Image3D = permute(FrFT1D_Uniform_Image3D, [3 2 1]);

for n = 1:N+1                                        % Traversing in the Z-axis
    Slice2D = FrFT1D_Uniform_Image3D(:,:,n);
    Polar2D = VectorizedCompute2DPolarDFT( Slice2D, M);
    CylindricalGrid(:,:,n) = Polar2D;
end
end

