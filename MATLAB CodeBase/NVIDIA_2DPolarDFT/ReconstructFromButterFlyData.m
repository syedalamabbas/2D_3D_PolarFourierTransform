function [ Image ] = ReconstructFromButterFlyData( FourierDataXaxisButterFlyGrid,  variableScalesY, uniformScaleX)
%RECONSTRUCTFROMBUTTERFLYDATA Summary of this function goes here
%   Detailed explanation goes here

[SizeX,~] = size(FourierDataXaxisButterFlyGrid);
Image = zeros(size(FourierDataXaxisButterFlyGrid));

%% order of operation is important , first we have to take variableScalesY only then we can work on uniform

for col = 1: SizeX
    Image(:, col) = VectorizedFrFT_Centered(FourierDataXaxisButterFlyGrid(:, col), -variableScalesY(col) ) ; % This is inverse FrFT
end

for row = 1: SizeX          
     Image(row,: ) = VectorizedFrFT_Centered(Image(row,: ).', -uniformScaleX ) ;      % This is inverse FrFT
end

% Image is now the final reconstruction
end

