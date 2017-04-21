function [ complex_point ] = FourierBasedTrigInterp( desired_point , poly_Coefficient , M ) % M = no of Angles
xCoordinate = desired_point(1);
yCoordinate = desired_point(2);
desired_theta = atan(yCoordinate/xCoordinate);
r = sqrt(xCoordinate^2+ yCoordinate^2);
spacing = [0:1:M-1, -M:1:-1]; % Set of angles
theta = spacing*desired_theta;
complex_point = sum (poly_Coefficient .* (r* exp(1i*theta)));
end

