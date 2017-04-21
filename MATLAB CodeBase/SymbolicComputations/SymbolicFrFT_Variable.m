function [ DFTsum ] = SymbolicFrFT_Variable( OneD_Vector_Symbols, beta, k , length_Vector, Size)
%SYMBOLICFFT_9ELEMENTS This does the computation of FrFT of 1 element
n = -(Size-1)/2:(Size-1)/2;
Map = exp (-1i* 2*pi* k * abs(k) * beta * n/length_Vector);
DFTsum = sum (OneD_Vector_Symbols .* Map);
end