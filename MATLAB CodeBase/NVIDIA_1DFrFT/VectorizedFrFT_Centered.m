function y=VectorizedFrFT_Centered(x,alpha)

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

[sizeX, ~] = size(x);

N =sizeX-1;   % N is always assumed even
n=[0:1:N, -N-1:1:-1]';  % Right

M = (0:1:N)';
k = (0:1:N)';
% Z = zeros(N+1,N+1);
Z = zeros(size(x));  
M = M - N/2; 


% double/triple padding ???
% n=[0:1:N, -N-1:1:-1, 0:1:N,  -N-1:1:-1]';
% Z = zeros(size([x;x;x]));

E_n=exp(-1i*pi*alpha*n.^2/(N+1)) ;          % Sequence as defined in the paper
% E_n(N+2:end) =Z;

% figure, 
% plot(real(E_n),'g')
% % hold on
% % plot(imag(E_n), 'r')
% % hold off
% title('Sequence plot')
PremultiplicationFactor=exp(1i*pi*M *N*alpha/(N+1));
PostmultiplicationFactor = exp(1i*pi*alpha*N*k/(N+1));


% Adding double for precision
% E_n= double(exp( double( -1i*pi*alpha*n.^2/(N+1))));          % Sequence as defined in the paper
% PremultiplicationFactor=double(exp( double( 1i*pi*M *N*alpha/(N+1))));
% PostmultiplicationFactor = double(exp( double(1i*pi*alpha*N*k/(N+1))));



x_tilde=  bsxfun(@times,x,PremultiplicationFactor);
x_tilde=[x_tilde; Z];

x_tilde= bsxfun(@times,x_tilde,E_n);

% These three operations constitute "fast convolution" or FFT based convolution
FirstFFT=fft(x_tilde);
SecondFFT=fft(conj(E_n));
interimY=ifft( bsxfun(@times,FirstFFT,SecondFFT));             % Computing Convolution here

interimY= bsxfun(@times,interimY,E_n);
y = bsxfun(@times,PostmultiplicationFactor,interimY(1:N+1,:,:));
return; 

