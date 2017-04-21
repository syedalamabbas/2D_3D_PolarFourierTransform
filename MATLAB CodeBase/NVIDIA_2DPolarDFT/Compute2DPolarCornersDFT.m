function [ PolarGridCorners,PolarGridCornersWeights, C ] = Compute2DPolarCornersDFT( inputImage,  noOfAngles )
% Compute2DPolarCornersDFT computes few additional measurement of corners that can be used along with COMPUTE2DRADONTRANSFORM 
% This function has been deprecated on 11/29/2015

I = inputImage;
[sizeX,~] =  size(I);
N = sizeX -1;      % N is always even
M = noOfAngles;    % M is also even

deltaTheta = 180/M;                     % Angular sampling rate
angles = deltaTheta:90-deltaTheta;      % Considering the first quadrant only  excluding the 0th and 90th for sure they dont and cant have any corner points

FirstQuadPoints = [];   
PolarGridCornersWeights = [];
UniformGridSpacing = -N/2:N/2;

for angle = angles
    if ( angle <= 45 )                   % BH lines
        newGridSpacing = UniformGridSpacing ./ cos(angle*pi/180);                     % BH   Since 0 < angle <= 45^o
    else
        newGridSpacing = UniformGridSpacing ./ sin(angle*pi/180);                     % BV   Since 45 < angle < 90^o
    end
    halfSpacing = 1:1:newGridSpacing(end);
    halfPoints = halfSpacing(N/2+1:end);
    
    if (~isempty(halfPoints))
        for point = halfPoints
            FirstQuadPoints = [FirstQuadPoints ; [point*cosd(angle),point*sind(angle)]];
            PolarGridCornersWeights = [PolarGridCornersWeights; sqrt(point/2)/(N+1)]  ;   % Analytical weights as given in Amir 2006 paper  sqrt(abs(-N/2:N/2)/2)/fullSize;
%             PolarGridCornersWeights = [PolarGridCornersWeights; pi/M* point/N^2]  ;             % Analytical weights as given in Markus 2007 Polar FFT Applied Computational Harmonic paper
        end
    end
end

%% Check the points if you have to
figure, 
hold on
scatter (FirstQuadPoints(:,1),FirstQuadPoints(:,2));
scatter (FirstQuadPoints(:,1),-FirstQuadPoints(:,2));
scatter (-FirstQuadPoints(:,1),FirstQuadPoints(:,2));
scatter (-FirstQuadPoints(:,1),-FirstQuadPoints(:,2));
hold off
axis equal; xlabel('x'); ylabel('y')

C = length(FirstQuadPoints);        % Number of corners  as computed in the first quadrant

% PolarGridCorners = zeros(1, C);          %  Polar grid corners: Four Quadrants (1st, 2nd, 3rd and 4th Rowwise) vs Corner Radial Data
PolarGridCorners = sym('asasd',[1, C]);    % Activate only for symbolic computations

for c =1:C
    PolarGridCorners(:,c) = CornerPoints_2DDFT( I, FirstQuadPoints(c,:) );
end
end

function [ FourDFT_Points ] = CornerPoints_2DDFT( inputImage, desiredPoint )
%   Computes the solution for a single point in the first quadrant  
%   directly but also computes the other 3 points in the other 3 quadrants
%   as a by product 
I = inputImage;
[sizeX,~] =  size(I);
N = sizeX -1;      % N is even

xIndex = desiredPoint(1);
yIndex = desiredPoint(2); 

indexes = -N/2:N/2;                                 % Range of indexes 
Map_x = exp(-1i*2*pi*xIndex*indexes/(N+1));         % Rows 
Map_y = exp(-1i*2*pi*yIndex*indexes' /(N+1));       % Columns

% lineData     = sum(bsxfun(@times,I,Map_y));             % One direct multiplication (N+1)^2 but very fast especially using GPU
lineData = sum(SymbolicMultiply(I,Map_y));                % Activate only for symbolic computations
ConjLineData = conj(lineData);

FirstQuadPoint  = sum (lineData .* Map_x);               % Need to check the directions
% SecondQuadPoint = conj(FirstQuadPoint);
% FourthQuadPoint = 0; sum (ConjLineData .* Map_x);
% ThirdQuadPoint  = conj(FourthQuadPoint);

% FourDFT_Points = [FirstQuadPoint; SecondQuadPoint; ThirdQuadPoint; FourthQuadPoint];
FourDFT_Points = [FirstQuadPoint]; % FourthQuadPoint];
end


function ColMultiplyImage = SymbolicMultiply(I,Map_y)
[sizeX,~] =  size(I);
N = sizeX -1;      % N is even
ColMultiplyImage = [];
for k = 1:N+1
    ColMultiplyImage = [ColMultiplyImage, I(:,k).* Map_y];
end
end