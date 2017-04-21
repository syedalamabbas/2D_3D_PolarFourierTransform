function [ ValueAtDesiredIndex ] = FrFTCenteredSingle( x, alpha, desiredIndex )
%FRFTCENTEREDSINGLE computes the value using the given index k F^{alpha}(k)
%using the data f(n)

[sizeX, sizeY] = size(x);

if (sizeX == 1)
    N = sizeY-1;  % N is even
    spacing = [-N/2:1:N/2];
else
    N = sizeX-1;  % N is even
    spacing = [-N/2:1:N/2]';
end

scaleMultiples = exp(-1i*2*pi*spacing *alpha* desiredIndex/ (N+1));
ValueAtDesiredIndex = sum (x .* scaleMultiples);
end

