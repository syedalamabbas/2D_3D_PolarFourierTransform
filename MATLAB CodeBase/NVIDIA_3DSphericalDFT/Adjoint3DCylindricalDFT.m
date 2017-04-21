function [ ImageAdjoint3D ] = Adjoint3DCylindricalDFT( CylindricalGrid )
%ADJOINT3DCYLINDRICALDFT  Transforms the fucntion using the prescribed
% scheme as discussed in the paper by Syed Alam Abbas
% Synopsis: Y=Adjoint3DCylindricalDFT(X)
%
% Inputs -  X  M x N+1 x (N+1) matrix in Cylindrical grid, (N is assumed to be even)           
% Outputs - Y  (N+1) x (N+1) x (N+1) matrix (full Cartesian grid)
%
% Example: 
% 
%   The following is a way to verify that this works as an adjoint -
%   choosing random X and Y, one must obtain that 
% The adjoint of a transformation A is a unique transformation A^H so that 
% <Ax,y> = <x, A^H y> for every x and y
% Adjoint operation
% N = 16;                                                % Always even
% M = N+2;                                               % Always even
% X = randn(N+1, N+1, N+1)+ 1i*randn(N+1, N+1, N+1);     % Complex matrix of Polar Grid
% Y = randn(M, N+1, N+1)  + 1i*randn(M, N+1, N+1);       % Complex matrix
% Ax = VectorizedCompute3DCylindricalDFT( X,  M );       % Output M x N+1 x N+1     Cylindrical grid
% AHy = Adjoint3DCylindricalDFT(Y);                      % Output N+1 x N+1 x N+1   Adjoint image
% % Checking Business
% prod1  = squeeze(dot(Ax,Y)) ;
% prod2  = squeeze(dot(X,AHy));
% fprintf('\n The result of  (<Ax,y> - <x, A^H y>) comparison is :')
% disp(abs(sum(sum(prod1)) -(sum(sum(prod2)))))
% Output is ---->>>>>   The result of  (<Ax,y> - <x, A^H y>) comparison is :   3.6265e-11
% Written on October 28th, 2015 by Syed Alam Abbas.
%====================================================================

% Z -coordinates vs  Polar slices: No of angles vs. Radial data
N =  size(CylindricalGrid,3) -1;
ImageAdjoint3D = zeros (N+1,N+1,N+1);                 % Output

for n = 1:N+1
    Polar2D = CylindricalGrid(:,:,n);
    Slice2D  = Adjoint2DPolarDFT( squeeze(Polar2D) );           % Using 2D adjoint
    ImageAdjoint3D(:,:,n) = Slice2D.';
end

ReorderedImage_OperateColumns = permute(ImageAdjoint3D, [3 2 1]);                     %  achieving with Swap Z-axis as columns now
FrFT1D_Uniform_Image3D = VectorizedFrFT_Centered(ReorderedImage_OperateColumns,-1);   %  Note the change in the sign
ImageAdjoint3D = permute(FrFT1D_Uniform_Image3D, [3 2 1]);                            %  Final Solution

end

