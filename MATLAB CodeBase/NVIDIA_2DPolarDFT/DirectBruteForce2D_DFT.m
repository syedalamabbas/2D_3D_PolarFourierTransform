function [ PolarGrid ] = DirectBruteForce2D_DFT( inputImage,  noOfAngles)
%DIRECTBRUTEFORCE2D_DFT Summary of this function goes here
%   Detailed explanation goes here
I = inputImage;
[sizeX,sizeY] =  size(I);
N = sizeX -1;      % N is even
M = noOfAngles;    % M is also even

PolarGrid = zeros(M, N+1);     %  Polar grid: No of angles vs. Radial data

deltaTheta = double(180/M);     % Angular sampling rate
angles = 0:deltaTheta:180-deltaTheta;

gridSpacing =  -N/2:N/2;
lineData = zeros(size(gridSpacing));

% figure,
% hold on
m =1;
for angle = angles
%     plot (gridSpacing*cosd(angle) ,gridSpacing*sind(angle), '-go')   % Debugging only
    for k =1:N+1
        lineData(k) = DirectSinglePoint_2DDFT( inputImage, [gridSpacing(k)*cosd(angle) ,gridSpacing(k)*sind(angle)  ] );
    end
    PolarGrid (m,:) = lineData;
    m = m+1;
end
% hold off  
end

 