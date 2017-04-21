function y=FrFT_Centered(x,alpha)

%=====================================================================
% Computes centered Fractional Fourier Transform according to the
% definition given in the paper: "An exact and fast computation of Discrete Fourier
% Transform for polar grid and spherical grid"
%
%           F^{\alpha}(k) = \sum\limits_{n=-N/2}^{N/2} f(n) e^{-j\frac{2\pi k\alpha n}{N+1}} ,
%           -N/2 \leq k \leq N/2
%
% For alpha=1 we get the regular FFT, and for alpha= -1 we get the regular
% IFFT up to a scale
%
% Synopsis: y=FRFT_Centered(x,alpha)   This is fully centered
%
% Inputs -    x - an N+1  - entry vector to be transformed
%                   alpha - the arbitrary scaling factor
%
% Outputs-  y - the transformed result as an N+1-entries vector
%=====================================================================

[sizeX, sizeY] = size(x);

if (sizeX == 1)
     N =sizeY-1;   % N is always assumed even
     n=[0:1:N, -N-1:1:-1];
     M = (0:1:N);
     k = (0:1:N);
     Z = zeros(1,N+1);
else
    N =sizeX-1;   % N is always assumed even
    n=[0:1:N, -N-1:1:-1]';
    M = (0:1:N)';
    k = (0:1:N)';
    Z = zeros(N+1,1);
end
M = M - N/2;


E_n=exp(-1i*pi*alpha*n.^2/(N+1)) ;          % Sequence as defined in the paper
PremultiplicationFactor=exp(1i*pi*M *N*alpha/(N+1));
PostmultiplicationFactor = exp(1i*pi*alpha*N*k/(N+1));

x_tilde=x.*PremultiplicationFactor;

if (sizeX == 1)
    x_tilde=[x_tilde, Z];
else
    x_tilde=[x_tilde; Z];
end

x_tilde=x_tilde.*E_n;


FirstFFT=fft(x_tilde);
SecondFFT=fft(conj(E_n));
interimY=ifft(FirstFFT.*SecondFFT);             % Computing Convolution here
interimY=interimY.*E_n;
y = PostmultiplicationFactor.*interimY(1:N+1);

return;

