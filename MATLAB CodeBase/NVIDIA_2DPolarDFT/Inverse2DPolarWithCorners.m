function [ ImageFinal ] = Inverse2DPolarWithCorners( Polar_Grid, PolarGridCorners, accuracy )
%INVERSE2DPOLARWITHCORNERS Summary of this function goes here
%   Detailed explanation goes here

if nargin==2,
    accuracy=1e-5;   % Default accuracy
end;

[M, sizeN_plus1 ] = size(Polar_Grid);
N = sizeN_plus1 -1;
ImageFinal = zeros(N+1,N+1);

fullSize  = sizeN_plus1;
W         = sqrt(abs(-N/2:N/2)/2)/fullSize;    % Fourier based preconditioner from Amir's paper
W(N/2+1)  =  sqrt(1/8)/fullSize;
W        = ones(M,1)*W ;


Delta=1;
count=0;
maxIterations = 5;
while Delta>accuracy && count < maxIterations,
    
    Err = W.* (Compute2DPolarDFT( ImageFinal,  M ) - Polar_Grid);                      % Residual from the polar system
    [CurrentPolarCorners,W2,~] = Compute2DPolarCornersDFT2( ImageFinal,  M );
    Err1 = W2 .* (CurrentPolarCorners - PolarGridCorners);                            % Residual from the polar corner system
    
    D1 = Adjoint2DPolarDFT(W.* Err).';
%     figure, imagesc(real(D1))
    D2 = Adjoint2DPolarCornersDFT2( W2 .* Err1, N,M ).';
%     figure, imagesc(real(D2))
    D =  D1 + D2;
    
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

end

