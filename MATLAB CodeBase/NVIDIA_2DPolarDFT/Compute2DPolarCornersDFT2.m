function [ PolarGridCorners,PolarGridCornersWeights, C ] = Compute2DPolarCornersDFT2( inputImage,  noOfAngles )
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
PolarGridCornersWeights = [PolarGridCornersWeights'; PolarGridCornersWeights'; PolarGridCornersWeights'; PolarGridCornersWeights'];  % Replicating for four corners

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
PolarGridCorners = zeros(4, C);     %  Polar grid corners: Four Quadrants (1st, 2nd, 3rd and 4th Rowwise) vs Corner Radial Data

for c =1:C
    point = FirstQuadPoints(c,:);
    xIndex = point(1);
    yIndex  = point(2);
    PolarGridCorners(:,c) = CornerPoints_2DDFT( I,  [xIndex , yIndex ] );
%     PolarGridCorners(1,c) = CornerPoints_2DDFT( I,  [xIndex , yIndex ] );  
%     PolarGridCorners(2,c) = CornerPoints_2DDFT( I,  [xIndex , -yIndex ] );  
%     PolarGridCorners(3,c) = CornerPoints_2DDFT( I,  [-xIndex , yIndex ] );  
%     PolarGridCorners(4,c) = CornerPoints_2DDFT( I,  [-xIndex , -yIndex ] );  
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

lineData     = sum(bsxfun(@times,I,Map_y));         % One direct multiplication (N+1)^2 but very fast especially using GPU
ConjLineData = conj(lineData);

%% Need to check the directions , naming/order doesnt matter the values are correct
FirstQuadPoint  = sum (lineData .* Map_x);           %  1 N+1 multiplication   
SecondQuadPoint = sum (ConjLineData .* Map_x);
ThirdQuadPoint  = conj(SecondQuadPoint);       %  1 N+1 multiplication   
FourthQuadPoint  = conj(FirstQuadPoint);

FourDFT_Points = [FirstQuadPoint; SecondQuadPoint; ThirdQuadPoint; FourthQuadPoint];
% FourDFT_Points = [FirstQuadPoint];
end


