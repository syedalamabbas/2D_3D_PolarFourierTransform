function [ matrix_Transform ] = CreateMatrixFrom3DGrid( charGrid, K,M,N, show)
%CREATEMATRIXFROMGRID creates the matrix transform operator
% N+1 x N+1 is the size of the image
% M is the number of angles

deltaTheta = 180/K;                    % Angular sampling rate theta
anglesTheta = 0:deltaTheta:180-deltaTheta;

deltaPhi = 180/M;                      % Angular sampling rate phi
anglesPhi = 0:deltaPhi:180-deltaPhi;


matrix_Transform = [];
gridSpacing = -N/2:N/2;
lineWidth = 1.7;
switch(charGrid)
%     case 'C'       % Pure Cartesian Grid
%         [CartesianGridX, CartesianGridY]= meshgrid(gridSpacing,gridSpacing);
%         for k = 1:length(gridSpacing)
%             for l = 1:length(gridSpacing)
%                 point = [CartesianGridX(k,l);CartesianGridY(k,l)];
%                 if(show)
%                         plot(point(1),point(2), '.r', 'LineWidth', lineWidth)
%                 end
%                 matrix_Transform = [ matrix_Transform; CreateRowFor3DTransform( N, point )];
%             end
%         end
    case 'P'         % Pure Polar Grid
        for k = 1:K
            angleTheta = anglesTheta(k);
            for m = 1:M
                anglePhi = anglesPhi (m);
                for n = 1:N+1
                    u = gridSpacing(n)*cosd(angleTheta)*cosd(anglePhi);
                    v = gridSpacing(n)*cosd(angleTheta)*sind(anglePhi);
                    w = gridSpacing(n)*sind(angleTheta);
                    desiredPoint = [ u, v, w ];
                    matrix_Transform = [ matrix_Transform; CreateRowFor3DTransform( N, desiredPoint )];
                end
            end
        end
%     case 'P_C'       % Polar grid with corners included
%         for angle = angles
%             if ( angle <= 45 || angle > 135 )
%                 newGridSpacing = gridSpacing ./ cos(angle*pi/180);       % BH
%                 if(angle > 135)
%                     halfSpacing = 1:1:newGridSpacing(1);
%                 else
%                     halfSpacing = 1:1:newGridSpacing(end);
%                 end
%             else if ( angle > 45 || angle <= 135 )
%                     newGridSpacing = gridSpacing ./ sin(angle*pi/180);    %BV
%                     halfSpacing = 1:1:newGridSpacing(end);
%                 end
%             end
%             newGridSpacing = [-fliplr(halfSpacing) 0 halfSpacing];
%             Line = [newGridSpacing*cosd(angle);newGridSpacing*sind(angle)];
%             if(show)
%                 plot(newGridSpacing*cosd(angle),newGridSpacing*sind(angle),'.r', 'LineWidth', lineWidth)
%             end
%             for k = 1:length(Line)
%                 point = Line(:,k);
%                 matrix_Transform = [ matrix_Transform; CreateRowFor3DTransform( N, point )];
%             end
%         end
%         
%     case 'E_P'         % Pure Polar Grid but with doubled Radial points
%         newGridSpacing = [-N/2:.5:N/2];
%         for angle = angles
%             Line = [newGridSpacing*cosd(angle);newGridSpacing*sind(angle)];
%             if(show)
%                 plot(newGridSpacing*cosd(angle),newGridSpacing*sind(angle), '.r', 'LineWidth', lineWidth)
%             end
%             for k = 1:length(Line)
%                 point = Line(:,k);
%                 matrix_Transform = [ matrix_Transform; CreateRowFor3DTransform( N, point )];
%             end
%         end
%     case 'E_PC'         % Polar Grid with corners included but with doubled Radial points
%         for angle = angles
%             if ( angle <= 45 || angle > 135 )
%                 newGridSpacing = gridSpacing ./ cos(angle*pi/180);       % BH
%                 if(angle > 135)
%                     halfSpacing = 1:1:newGridSpacing(1);
%                 else
%                     halfSpacing = 1:1:newGridSpacing(end);
%                 end
%             else if ( angle > 45 || angle <= 135 )
%                     newGridSpacing = gridSpacing ./ sin(angle*pi/180);    %BV
%                     halfSpacing = 1:1:newGridSpacing(end);
%                 end
%             end
%             newGridSpacing = [-halfSpacing(end):.5:halfSpacing(end) ];
%             Line = [newGridSpacing*cosd(angle);newGridSpacing*sind(angle)];
%             if(show)
%                 plot(newGridSpacing*cosd(angle),newGridSpacing*sind(angle),'.r', 'LineWidth', lineWidth)
%             end
%             
%             for k = 1:length(Line)
%                 point = Line(:,k);
%                 matrix_Transform = [ matrix_Transform; CreateRowFor3DTransform( N, point )];
%             end
%         end
%         
%     case 'P_CartC'         % Pure Polar Grid but with doubled Radial points
%         newGridSpacing = [-N/2:1:N/2];
%         PolarGridX = [];
%         PolarGridY = [];
%         for angle = angles
%             Line = [newGridSpacing*cosd(angle);newGridSpacing*sind(angle)];
%             if(show)
%                 plot(newGridSpacing*cosd(angle),newGridSpacing*sind(angle),'.r' , 'LineWidth', lineWidth)
%             end
%             for k = 1:length(Line)
%                 point = Line(:,k);
%                 matrix_Transform = [ matrix_Transform; CreateRowFor3DTransform( N, point )];
%             end
%         end
%         
%          newGridSpacing = [-N/2:1:N/2];
%         [CartesianGridX, CartesianGridY]= meshgrid(newGridSpacing,newGridSpacing);
%         matrix_Transform2 = [];
%         for k = 1:length(newGridSpacing)
%             for l = 1:length(newGridSpacing)
%                 point = [CartesianGridX(k,l);CartesianGridY(k,l)];
%                 matrix_Transform2 = [ matrix_Transform2; CreateRowFor3DTransform( N, point )];
%             end
%         end
%         Weight=(CartesianGridX.^2 + CartesianGridY.^2 > (N/2)^2); % All points outside the circle radius N/2 are weight 1 ,rest are zero
% %         figure, imagesc(Weight)
%         
%         if(show)
%             for k = 1:length(newGridSpacing)
%                 for l = 1:length(newGridSpacing)
%                     if (Weight(k,l))
%                         point = [CartesianGridX(k,l);CartesianGridY(k,l)];
%                         plot(point(1),point(2), '+b', 'LineWidth', lineWidth)
%                     end
%                 end
%             end
%         end
%                 
%         Pos=find(Weight(:));
%         matrix_Transform = [matrix_Transform; matrix_Transform2(Pos,:) ];
%         
%         case 'EP_CartC'         % Polar Grid but with doubled Radial points and Cartesian corners ! Best configuration
%             factor = .5;          
%         newGridSpacing = [-N/2:factor:N/2];
%         PolarGridX = [];
%         PolarGridY = [];
%         for angle = angles
%             Line = [newGridSpacing*cosd(angle);newGridSpacing*sind(angle)];
%             if(show)
%                 plot(newGridSpacing*cosd(angle),newGridSpacing*sind(angle),'.r', 'LineWidth', lineWidth )
%             end
%             for k = 1:length(Line)
%                 point = Line(:,k);
%                 matrix_Transform = [ matrix_Transform; CreateRowFor3DTransform( N, point )];
%             end
%         end
%         
%         [CartesianGridX, CartesianGridY]= meshgrid(newGridSpacing,newGridSpacing);
%         matrix_Transform2 = [];
%         for k = 1:length(newGridSpacing)
%             for l = 1:length(newGridSpacing)
%                 point = [CartesianGridX(k,l);CartesianGridY(k,l)];
%                 matrix_Transform2 = [ matrix_Transform2; CreateRowFor3DTransform( N, point )];
%             end
%         end
%         Weight=(CartesianGridX.^2 + CartesianGridY.^2 > (N/2)^2); % All points outside the circle radius N/2 are weight 1 ,rest are zero
% %         figure, imagesc(Weight)
%         
%         if(show)
%             for k = 1:length(newGridSpacing)
%                 for l = 1:length(newGridSpacing)
%                     if (Weight(k,l))
%                         point = [CartesianGridX(k,l);CartesianGridY(k,l)];
%                         plot(point(1),point(2), '+b', 'LineWidth', lineWidth)
%                     end
%                 end
%             end
%         end
%                 
%         Pos=find(Weight(:));
%         matrix_Transform = [matrix_Transform; matrix_Transform2(Pos,:) ];
%         
end

end


function [ Multiplication_Row ] = CreateRowFor3DTransform( N, desiredPoint )
%   Computes the solution for a single point in the first quadrant

xIndex = desiredPoint(1);             % The desired point is in the grid
yIndex = desiredPoint(2);
zIndex = desiredPoint(3);

indexes = -N/2:N/2;                                 % Range of indexes
Map_x = exp(-1i*2*pi*xIndex*indexes/(N+1));         % Rows
Map_y = exp(-1i*2*pi*yIndex*indexes' /(N+1));       % Columns
Map_z = exp(-1i*2*pi*zIndex*indexes' /(N+1));       % depth

[GridX, GridY, GridZ] = meshgrid(Map_x,Map_y, Map_z);

Grid = GridX .* GridY .* GridZ;
Multiplication_Row = reshape(Grid, [1, (N+1)^3]);
end

