function [ DFT_Point ] = DirectSinglePoint_2DDFT( inputImage, desiredPoint )
%COMPUTEFORSINGLEPOINT Summary of this function goes here
%   Computing the solution for a single point using row and column
%   differently scaled FrFT
I = inputImage;
[sizeX,sizeY] =  size(I);
N = sizeX -1;      % N is even

xIndex = desiredPoint(1);
yIndex = desiredPoint(2);

lineData = zeros(1, N+1);
y = -N/2: N/2;
Map_y = exp(-1i*2*pi*yIndex*y/(N+1));
Map_x = exp(-1i*2*pi*xIndex*y/(N+1));
                                             % Fully direct computations
for x = 1: N+1     % x
    col = (inputImage(:,x));
    lineData(x) = sum( col' .* Map_y);
end
DFT_Point =  sum (lineData .* Map_x);
end
 
