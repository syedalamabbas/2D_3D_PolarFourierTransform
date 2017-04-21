clc
clear
close all

addpath('..\NVIDIA_2DPolarDFT')
%% Brute force direct computation of 3-level BTTB matrix 
for N = 18; %[14:2:18]
    K = N+2;
    M = N + 2;
    HugeTransformMatrix = CreateMatrixFrom3DGrid( 'P', K,M,N );
end

BTTB_3Level_FullMat = HugeTransformMatrix'*HugeTransformMatrix;

figure,
imagesc(real(BTTB_3Level_FullMat))
title('What kind of Toeplitz is this ?- 3Level BTTB ')

%% Verify just 2D multiplication of BTTB matrix-vector product
BTTBColumnsIn2DShape = reshape(BTTB_3Level_FullMat(1:(N+1)^2,1), [N+1,N+1]);
figure, imagesc(real(BTTBColumnsIn2DShape))
title('Single column BTTB matrix arranged in 2D shape')
figure, imagesc(real(BTTB_3Level_FullMat(1:(N+1)^2,1:(N+1)^2)));
title('Full BTTB inner matrix')

x_BTTB_1D = rand((N+1)^2,1);                                                % Random input 
b_directMul_BTTB  = BTTB_3Level_FullMat(1:(N+1)^2,1:(N+1)^2)*x_BTTB_1D;     % Direct multiplication    

BCCB_EigenValues = fft2( makeBCCBFromBTTB_2Dextension(BTTBColumnsIn2DShape,N));
b_iterate_Inner_1D = BTTB_2D_Matrix_Vector(x_BTTB_1D,BCCB_EigenValues,N );      % Special structured Matrix multiplication

figure, plot (1:(N+1)^2, real(b_directMul_BTTB), 1:(N+1)^2, real(b_iterate_Inner_1D),'+' )   % Comparison
xlabel('signal, N \rightarrow')
ylabel('Intensity')
legend('Direct O(N^4)','With Structured Matrix O(N^2 log(N))');
title('BTTB 2 Level comparison') 
grid on
axis tight
disp(max(real(b_directMul_BTTB - b_iterate_Inner_1D)))
%% Verification of MAtrix product multiplication Ax = b , where A  is a 3-level BTTB matrix computed as above
x_iterate_1D = rand((N+1)^3,1);
tic 
b_directMul_BTTB_3D = BTTB_3Level_FullMat * x_iterate_1D;
toc 

BTTB_ColumnsIn3DShape = reshape(BTTB_3Level_FullMat(1,:), [N+1,N+1,N+1]);
BCCB_3D_EigenValues = zeros(2*N,2*N,N+1);
for k = 1:N+1
    BCCB_3D_EigenValues(:,:,k) = fft2(makeBCCBFromBTTB_2Dextension(BTTB_ColumnsIn3DShape(:,:,k),N));  % For a BCCB matrix  we have the spectral decomposition, A = F* Eig F, F* is the ifft2, F is fft2
end
tic 
b_iterate_1D =  BTTB_3D_Matrix_Vector(x_iterate_1D,BCCB_3D_EigenValues,N );
toc
figure, plot (1:(N+1)^3, real(b_directMul_BTTB_3D), 1:(N+1)^3, real(b_iterate_1D),'+' )   % Comparison
xlabel('signal, N \rightarrow')
ylabel('Intensity')
legend('Direct  O(N^6)','With Structured Matrix, O(N^3 log(N))');
title('BTTB 3 Level comparison')
grid on
axis tight

disp(max(real(b_directMul_BTTB_3D - b_iterate_1D)))