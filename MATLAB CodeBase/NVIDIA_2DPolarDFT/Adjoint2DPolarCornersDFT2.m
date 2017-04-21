function [ ImageAdjoint ] = Adjoint2DPolarCornersDFT2( PolarGridCorners, N,M )
%ADJOINT2DPOLARCORNERSDFT  computes the adjoint or the reverse of the
%Compute2DPolarCornersDFT operation
% This function has been deprecated on 11/29/2015

deltaTheta = 180/M;                     % Angular sampling rate
angles = deltaTheta:90-deltaTheta;      % Considering the first quadrant only  excluding the 0th and 90th for sure they dont and cant have any corner points

FirstQuadPoints = [];
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
        end
    end
end

C = length(FirstQuadPoints);        % Number of corners  as computed in the first quadrant
ImageAdjoint = zeros(N+1);
for c =1:C
    point = FirstQuadPoints(c,:);
    xIndex = point(1);
    yIndex  = point(2);
    ImageAdjoint = ImageAdjoint + CornerPoints_2DAdjointDFT(PolarGridCorners(1,c),N, [xIndex , yIndex ] ) ;       % Collecting contributions from all corner points
    ImageAdjoint = ImageAdjoint + CornerPoints_2DAdjointDFT(PolarGridCorners(2,c),N, [xIndex , -yIndex ] ) ; 
    ImageAdjoint = ImageAdjoint + CornerPoints_2DAdjointDFT(PolarGridCorners(3,c),N, [-xIndex , yIndex ] ) ; 
    ImageAdjoint = ImageAdjoint + CornerPoints_2DAdjointDFT(PolarGridCorners(4,c),N, [-xIndex , -yIndex ] ) ; 
end
end


function [ Image_Contribution ] = CornerPoints_2DAdjointDFT( PolarGridCornersFourPoints,N, desiredPoint )

xIndex = desiredPoint(1);
yIndex = desiredPoint(2);

indexes = -N/2:N/2;                                 % Range of indexes
Map_x = exp(+1i*2*pi*xIndex*indexes' /(N+1));         % Rows  Observe the sign change from -ve to +ve
Map_y = exp(+1i*2*pi*yIndex*indexes /(N+1));       % Columns  Observe the sign change from -ve to +ve


%% Reversing operation      -------> FirstQuadPoint  = sum (lineData .* Map_x);
FirstQuadPoint  = PolarGridCornersFourPoints;
%   SecondQuadPoint = PolarGridCornersFourPoints(2);
%   ThirdQuadPoint  = PolarGridCornersFourPoints(3);
%   FourthQuadPoint = PolarGridCornersFourPoints(2);

lineData = FirstQuadPoint .* Map_y;            % New Row formed from a the corner points
% lineData2 = SecondQuadPoint .* conj(Map_y);
%  lineData3 = ThirdQuadPoint .* Map_y;
%  lineData4 = FourthQuadPoint .* conj(Map_y);

image1 = bsxfun(@times,ones(N+1,1) *lineData, Map_x);
% image2 = bsxfun(@times,ones(N+1,1) *lineData2, Map_x);
% image3 = bsxfun(@times,ones(N+1,1) *lineData3, conj(Map_x));
% image4 = bsxfun(@times,ones(N+1,1) *lineData4, conj(Map_x));

%% Reversing operation     ------->     lineData     = sum(bsxfun(@times,I,Map_y));
% Image_Contribution = image1+ image2+ image3 + image4; % Replicating the rows to form an (N+1) x (N+1) image
Image_Contribution = image1;
end
