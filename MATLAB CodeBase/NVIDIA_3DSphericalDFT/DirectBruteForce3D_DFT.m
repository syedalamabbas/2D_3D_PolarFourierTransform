function [ SphericalGrid ] = DirectBruteForce3D_DFT( inputVolume,  noOfAnglesTheta, noOfAnglesPhi)
%DIRECTBRUTEFORCE2D_DFT Summary of this function goes here

I = inputVolume;                        % 3D cube image
[sizeX, ~, ~] =  size(I);
N = sizeX -1;                          % N is even
K = noOfAnglesTheta;                   % K is also even it is for theta
M = noOfAnglesPhi;                     % M is also even it is for phi

SphericalGrid = zeros( K,M, N+1);      % Angle phi vs  Polar slices: No of angles vs. Radial data

deltaTheta = 180/K;                    % Angular sampling rate theta
anglesTheta = 0:deltaTheta:180-deltaTheta;

deltaPhi = 180/M;                      % Angular sampling rate phi
anglesPhi = 0:deltaPhi:180-deltaPhi;

gridSpacing =  -N/2:N/2;
debug = 0;
for k = 1:K
        angleTheta = anglesTheta(k);
    for m = 1:M
         anglePhi = anglesPhi (m);
        for n = 1:N+1
            
%             desiredPoint = [gridSpacing(n)*sind(angleTheta)*cosd(anglePhi),gridSpacing(n)*sind(angleTheta)*sind(anglePhi),gridSpacing(n)*cosd(angleTheta)];
            u = gridSpacing(n)*cosd(angleTheta)*cosd(anglePhi);
            v = gridSpacing(n)*cosd(angleTheta)*sind(anglePhi);
            w = gridSpacing(n)*sind(angleTheta);
            desiredPoint = [ u, v, w ];
%             if (1 == CheckPairsXX(k,m,K,M))   %% This condition is for checking only XX block
%             if(1 == CheckPairsXZ(k,m,K,M))
                SphericalGrid(k,m,n)= DirectSinglePoint_3DDFT( I, desiredPoint );
%             end
        end   
        if (debug && k == 2 && m == 2) 
            figure, scatter3(gridSpacing*cosd(angleTheta)*cosd(anglePhi),gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta), 'filled');
            xlabel('x')
            ylabel('y')
            zlabel('z')
            hold on
            plot3(gridSpacing*cosd(angleTheta)*cosd(anglePhi),gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta), 'g','LineWidth', 1.7)
            hold off 
        end 
    end 
end 
end

