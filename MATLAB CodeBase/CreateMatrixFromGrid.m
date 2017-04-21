function [ matrix_Transform ] = CreateMatrixFromGrid( charGrid, N ,M , show)
%CREATEMATRIXFROMGRID creates the matrix transform operator
% N+1 x N+1 is the size of the image
% M is the number of angles

deltaTheta = 180/M;     % Angular sampling rate
angles = [0:deltaTheta:180-deltaTheta];
matrix_Transform = [];
gridSpacing = [-N/2:N/2];
lineWidth = 1.7;
switch(charGrid)
    case 'C'       % Pure Cartesian Grid
        [CartesianGridX, CartesianGridY]= meshgrid(gridSpacing,gridSpacing);
        for k = 1:length(gridSpacing)
            for l = 1:length(gridSpacing)
                point = [CartesianGridX(k,l);CartesianGridY(k,l)];
                if(show)
                        plot(point(1),point(2), '.r', 'LineWidth', lineWidth)
                end
                matrix_Transform = [ matrix_Transform; CreateRowForTransform( N, point )];
            end
        end
    case 'P'         % Pure Polar Grid
        for angle = angles
            Line = [gridSpacing*cosd(angle);gridSpacing*sind(angle)];
            if(show)
                plot(gridSpacing*cosd(angle),gridSpacing*sind(angle), '.r', 'LineWidth', lineWidth)
            end
            for k = 1:length(Line)
                point = Line(:,k);
                matrix_Transform = [ matrix_Transform; CreateRowForTransform( N, point )];
            end
        end
    case 'P_C'       % Polar grid with corners included
        for angle = angles
            if ( angle <= 45 || angle > 135 )
                newGridSpacing = gridSpacing ./ cos(angle*pi/180);       % BH
                if(angle > 135)
                    halfSpacing = 1:1:newGridSpacing(1);
                else
                    halfSpacing = 1:1:newGridSpacing(end);
                end
            else if ( angle > 45 || angle <= 135 )
                    newGridSpacing = gridSpacing ./ sin(angle*pi/180);    %BV
                    halfSpacing = 1:1:newGridSpacing(end);
                end
            end
            newGridSpacing = [-fliplr(halfSpacing) 0 halfSpacing];
            Line = [newGridSpacing*cosd(angle);newGridSpacing*sind(angle)];
            if(show)
                plot(newGridSpacing*cosd(angle),newGridSpacing*sind(angle),'.r', 'LineWidth', lineWidth)
            end
            for k = 1:length(Line)
                point = Line(:,k);
                matrix_Transform = [ matrix_Transform; CreateRowForTransform( N, point )];
            end
        end
        
    case 'E_P'         % Pure Polar Grid but with doubled Radial points
        newGridSpacing = [-N/2:.5:N/2];
        for angle = angles
            Line = [newGridSpacing*cosd(angle);newGridSpacing*sind(angle)];
            if(show)
                plot(newGridSpacing*cosd(angle),newGridSpacing*sind(angle), '.r', 'LineWidth', lineWidth)
            end
            for k = 1:length(Line)
                point = Line(:,k);
                matrix_Transform = [ matrix_Transform; CreateRowForTransform( N, point )];
            end
        end
    case 'E_PC'         % Polar Grid with corners included but with doubled Radial points
        for angle = angles
            if ( angle <= 45 || angle > 135 )
                newGridSpacing = gridSpacing ./ cos(angle*pi/180);       % BH
                if(angle > 135)
                    halfSpacing = 1:1:newGridSpacing(1);
                else
                    halfSpacing = 1:1:newGridSpacing(end);
                end
            else if ( angle > 45 || angle <= 135 )
                    newGridSpacing = gridSpacing ./ sin(angle*pi/180);    %BV
                    halfSpacing = 1:1:newGridSpacing(end);
                end
            end
            newGridSpacing = [-halfSpacing(end):.5:halfSpacing(end) ];
            Line = [newGridSpacing*cosd(angle);newGridSpacing*sind(angle)];
            if(show)
                plot(newGridSpacing*cosd(angle),newGridSpacing*sind(angle),'.r', 'LineWidth', lineWidth)
            end
            
            for k = 1:length(Line)
                point = Line(:,k);
                matrix_Transform = [ matrix_Transform; CreateRowForTransform( N, point )];
            end
        end
        
    case 'P_CartC'         % Pure Polar Grid but with doubled Radial points
        newGridSpacing = [-N/2:1:N/2];
        PolarGridX = [];
        PolarGridY = [];
        for angle = angles
            Line = [newGridSpacing*cosd(angle);newGridSpacing*sind(angle)];
            if(show)
                plot(newGridSpacing*cosd(angle),newGridSpacing*sind(angle),'.r' , 'LineWidth', lineWidth)
            end
            for k = 1:length(Line)
                point = Line(:,k);
                matrix_Transform = [ matrix_Transform; CreateRowForTransform( N, point )];
            end
        end
        
         newGridSpacing = [-N/2:1:N/2];
        [CartesianGridX, CartesianGridY]= meshgrid(newGridSpacing,newGridSpacing);
        matrix_Transform2 = [];
        for k = 1:length(newGridSpacing)
            for l = 1:length(newGridSpacing)
                point = [CartesianGridX(k,l);CartesianGridY(k,l)];
                matrix_Transform2 = [ matrix_Transform2; CreateRowForTransform( N, point )];
            end
        end
        Weight=(CartesianGridX.^2 + CartesianGridY.^2 > (N/2)^2); % All points outside the circle radius N/2 are weight 1 ,rest are zero
%         figure, imagesc(Weight)
        
        if(show)
            for k = 1:length(newGridSpacing)
                for l = 1:length(newGridSpacing)
                    if (Weight(k,l))
                        point = [CartesianGridX(k,l);CartesianGridY(k,l)];
                        plot(point(1),point(2), '+b', 'LineWidth', lineWidth)
                    end
                end
            end
        end
                
        Pos=find(Weight(:));
        matrix_Transform = [matrix_Transform; matrix_Transform2(Pos,:) ];
        
        case 'EP_CartC'         % Polar Grid but with doubled Radial points and Cartesian corners ! Best configuration
            factor = .5;          
        newGridSpacing = [-N/2:factor:N/2];
        PolarGridX = [];
        PolarGridY = [];
        for angle = angles
            Line = [newGridSpacing*cosd(angle);newGridSpacing*sind(angle)];
            if(show)
                plot(newGridSpacing*cosd(angle),newGridSpacing*sind(angle),'.r', 'LineWidth', lineWidth )
            end
            for k = 1:length(Line)
                point = Line(:,k);
                matrix_Transform = [ matrix_Transform; CreateRowForTransform( N, point )];
            end
        end
        
        [CartesianGridX, CartesianGridY]= meshgrid(newGridSpacing,newGridSpacing);
        matrix_Transform2 = [];
        for k = 1:length(newGridSpacing)
            for l = 1:length(newGridSpacing)
                point = [CartesianGridX(k,l);CartesianGridY(k,l)];
                matrix_Transform2 = [ matrix_Transform2; CreateRowForTransform( N, point )];
            end
        end
        Weight=(CartesianGridX.^2 + CartesianGridY.^2 > (N/2)^2); % All points outside the circle radius N/2 are weight 1 ,rest are zero
%         figure, imagesc(Weight)
        
        if(show)
            for k = 1:length(newGridSpacing)
                for l = 1:length(newGridSpacing)
                    if (Weight(k,l))
                        point = [CartesianGridX(k,l);CartesianGridY(k,l)];
                        plot(point(1),point(2), '+b', 'LineWidth', lineWidth)
                    end
                end
            end
        end
                
        Pos=find(Weight(:));
        matrix_Transform = [matrix_Transform; matrix_Transform2(Pos,:) ];
        
end

end


function [ Multiplication_Row ] = CreateRowForTransform( N, desiredPoint )
%   Computes the solution for a single point in the first quadrant

xIndex = desiredPoint(1);             % The desired point is in the grid
yIndex = desiredPoint(2);

indexes = -N/2:N/2;                                 % Range of indexes
Map_x = exp(-1i*2*pi*xIndex*indexes/(N+1));         % Rows
Map_y = exp(-1i*2*pi*yIndex*indexes' /(N+1));       % Columns

[GridX, GridY] = meshgrid(Map_x,Map_y);
% Grid1 = Map_y * ones(1, N+1);
% Grid2 = ones( N+1,1) * Map_x;

Grid = GridX .* GridY;
Multiplication_Row = reshape(Grid, [1, (N+1)^2]);
end

