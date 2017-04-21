function y= GPUVectorizedFrFT_Centered(I, Z, E_n, PremultiplicationFactor, PostmultiplicationFactor )

%=====================================================================
% Computes centered Fractional Fourier Transform according to the
% definition given in the paper: "An exact and fast computation of Discrete Fourier
% Transform for polar grid"
%
%           F^{\alpha}(k) = \sum\limits_{n=-N/2}^{N/2} f(n) e^{-j\frac{2\pi k\alpha n}{N+1}} ,
%           -N/2 \leq k \leq N/2
%
% For alpha=1 we get the regular FFT, and for alpha= -1 we get the regular
% IFFT up to a scale
%
% Synopsis: y=VectorizedFrFT_Centered(x,alpha)   This is fully centered
%
% Inputs -  I - an N+1 x N+1  - entry vector to be transformed
%                   alpha - the arbitrary scaling factor
%
% Outputs-  y - the transformed result as an N+1-entries vector
%=====================================================================
[sizeX,~] =  size(I);
N = sizeX -1;      % N is even 

I_tilde=  bsxfun(@times,I,PremultiplicationFactor);
I_tilde=[I_tilde; Z];

I_tilde= bsxfun(@times,I_tilde,E_n);
FirstFFT=fft(I_tilde);
SecondFFT=fft(conj(E_n));
interimY=ifft( bsxfun(@times,FirstFFT,SecondFFT));             % Computing Convolution here
interimY= bsxfun(@times,interimY,E_n);
y = bsxfun(@times,interimY(1:N+1,:,:),PostmultiplicationFactor);
return; 
 
