function [ outputVolume ] = DirectBruteForce3D_AdjointDFT( SphericalGrid )
%DIRECTBRUTEFORCE3D_ADJOINTDFT Summary of this function goes here
%   Detailed explanation goes here
[~,~,SizeX] = size(SphericalGrid);                        % 3D spherical grid
N = SizeX -1;
gridSpacing =  -N/2:N/2;

outputVolume = zeros(N+1,N+1,N+1);                    % Assume its a cube
for r = 1:N+1
    for c = 1:N+1
        for d = 1:N+1
            desiredPoint = [ gridSpacing(r) gridSpacing(c) gridSpacing(d) ];
            outputVolume(r,c,d)= DirectSinglePoint_3DAdjointDFT( SphericalGrid, desiredPoint );
        end
    end
end
end

