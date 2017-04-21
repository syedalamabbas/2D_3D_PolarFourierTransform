function [ DFTsum ] = SymbolicFrFT_Uniform( OneD_Vector_Symbols, alpha, k ,  Size)
%SYMBOLICFFT_9ELEMENTS This does the computation of FrFT of 1 element
n = -(Size-1)/2:(Size-1)/2;
Map = exp (-1i* 2*pi* k * alpha* n/Size);
DFTsum = sum (OneD_Vector_Symbols .* Map);
end

