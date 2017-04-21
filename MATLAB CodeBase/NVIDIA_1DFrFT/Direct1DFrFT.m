function [ y ] = Direct1DFrFT( x,alpha )
%DIRECT1DFRFTCOMPUTE Summary of this function goes here
%   Detailed explanation goes here
N=length(x)-1;   % N is even
n = [-N/2:N/2];
y = zeros(size(x));
m = 1;
for k =-N/2:N/2             % Direct multiplication
    E_n = exp(-2i*pi* alpha* k * n / (N+1));
    y(m) = sum(x .* E_n);
    m = m+1;
end
end