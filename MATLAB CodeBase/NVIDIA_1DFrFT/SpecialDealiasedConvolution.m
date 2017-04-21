function [ cconv ] = SpecialDealiasedConvolution( f, g )
%SPECIALDEALIASEDCONVOLUTION This code does a highly specialized inplace
%convolution operation using ffts, also known as fast convolution
% Due to the implicit nature of the solution
%====================================================================
% Written on May 10th, 2016 by Syed Alam Abbas.
%====================================================================

% Input vector f - Complex
% Input vector g - Complex
% Output vector cconv - Complex

m = length(f);

u = fft(f);
v = fft(g);
u = bsxfun(@times,u,v);

mul_factor = exp(-2*pi*1i*(0:1:m-1)'/(2*m));
f = bsxfun(@times,f,mul_factor );
g = bsxfun(@times,g,mul_factor );

v = fft(f);
f = fft(g);
v = bsxfun(@times,v,f);

f = ifft(u);
u = ifft(v);

temp = bsxfun(@times,u,conj(mul_factor) );
f = f + temp;

cconv = f/(2);            % Final answer

end

