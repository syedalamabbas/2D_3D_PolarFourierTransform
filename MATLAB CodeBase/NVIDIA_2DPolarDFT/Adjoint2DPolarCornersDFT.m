function [ ImageAdjoint ] = Adjoint2DPolarCornersDFT( PolarGridCorners, N,M )
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
    ImageAdjoint = ImageAdjoint + CornerPoints_2DAdjointDFT(PolarGridCorners(:,c),N, FirstQuadPoints(c,:) ) ;       % Collecting contributions from all corner points
end
end


function [ Image_Contribution ] = CornerPoints_2DAdjointDFT( PolarGridCornersFourPoints,N, desiredPoint )

xIndex = desiredPoint(1);
yIndex = desiredPoint(2);

indexes = -N/2:N/2;                                 % Range of indexes
Map_x = exp(+1i*2*pi*xIndex*indexes/(N+1));         % Rows  Observe the sign change from -ve to +ve
Map_y = exp(+1i*2*pi*yIndex*indexes' /(N+1));       % Columns  Observe the sign change from -ve to +ve


%% Reversing operation      -------> FirstQuadPoint  = sum (lineData .* Map_x);
FirstQuadPoint  = PolarGridCornersFourPoints(1);
% SecondQuadPoint = PolarGridCornersFourPoints(2);
% ThirdQuadPoint  = PolarGridCornersFourPoints(3);
% FourthQuadPoint = PolarGridCornersFourPoints(2);

lineData = FirstQuadPoint .* Map_x;            % New Row formed from a the corner points
% lineData = lineData + (conj(FirstQuadPoint)*ones(1,N+1).* conj(Map_x));
% ConjlineData = FourthQuadPoint*ones(1,N+1);   % These are conjugate points
% ConjlineData = ConjlineData .* ( Map_x);            % New Row formed from a the corner points
% ConjlineData = ConjlineData + ((ThirdQuadPoint)*ones(1,N+1) .* conj( Map_x));       

%% Reversing operation     ------->     lineData     = sum(bsxfun(@times,I,Map_y));
% Image_Contribution =  bsxfun(@times,ones(N+1,1) * lineData,Map_y); % Replicating the rows to form an (N+1) x (N+1) image
% Image_Contribution  = Image_Contribution  + bsxfun(@times,ones(N+1,1) * ConjlineData ,conj(Map_y));

%% Activate only for symbolic computations 
Image_Contribution =  SymbolicMultiply( ones(N+1,1)* lineData  ,(Map_y)) ;
% Image_Contribution =  SymbolicMultiply( lineData.' * ones(1,N+1) ,(Map_y)) ;
% Image_Contribution  = Image_Contribution  +  SymbolicMultiply(ones(N+1,1) * ConjlineData,(Map_y));   % Swapped
end


function ColMultiplyImage = SymbolicMultiply(I,Map_y)        % Use this only for symbolic computations 
[sizeX,~] =  size(I);
N = sizeX -1;      % N is even
ColMultiplyImage = [];
for k = 1:N+1
    ColMultiplyImage = [ColMultiplyImage, I(:,k).* Map_y];
end
end
