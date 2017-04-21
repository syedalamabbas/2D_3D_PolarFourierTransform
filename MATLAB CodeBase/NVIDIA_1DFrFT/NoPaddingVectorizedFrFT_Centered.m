function y=NoPaddingVectorizedFrFT_Centered(x,alpha)
% This function has a special property that it has no padding !
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
% Inputs -  x - an N+1  - entry vector or a column  [ (N+1) x 1 ] to be transformed
%                   alpha - the arbitrary scaling factor
%
% Outputs-  y - the transformed result as an N+1-entries vector
%=====================================================================

%====================================================================
% Written on May 10th, 2016 by Syed Alam Abbas.
%====================================================================


[sizeX, ~] = size(x);
N =sizeX-1;   % N is always assumed even
n=[0:1:N, -N-1:1:-1]';  % Right
% n=(0:1:N)';  % test
M = (0:1:N)';
k = (0:1:N)';
M = M - N/2; 
Z = zeros(size(x));  

%% Can be precomputed 
E_n=exp(-1i*pi*alpha*n.^2/(N+1)) ;          % Sequence as defined in the paper
PremultiplicationFactor=exp(1i*pi*M *N*alpha/(N+1));
PostmultiplicationFactor = exp(1i*pi*alpha*N*k/(N+1));

%% Convolution main and some multiplication factors
x_tilde=  bsxfun(@times,x,PremultiplicationFactor);
x_tilde=[x_tilde; Z];
x_tilde= bsxfun(@times,x_tilde,E_n);
interimY = SpecialDealiasedConvolution( x_tilde, conj(E_n) );   % Computing Convolution here
interimY= bsxfun(@times,interimY,E_n);
y = bsxfun(@times,PostmultiplicationFactor,interimY(1:N+1,:,:));

return; 

