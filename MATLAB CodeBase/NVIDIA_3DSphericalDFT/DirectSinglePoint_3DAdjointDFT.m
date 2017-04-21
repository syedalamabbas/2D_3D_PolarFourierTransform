function [ AdjointDFT_Point ] = DirectSinglePoint_3DAdjointDFT( SphericalGrid, desiredPoint )
%DIRECTSINGLEPOINT_3DADJOINTDFT Summary of this function goes here
%   Detailed explanation goes here

[K,M,sizeX] =  size(SphericalGrid);             % Spherical grid K x M x (N+1) with 

N = sizeX -1;      % N is even

r = desiredPoint(1);
c = desiredPoint(2);
d = desiredPoint(3);

deltaTheta = 180/K;                    % Angular sampling rate theta
anglesTheta = 0:deltaTheta:180-deltaTheta;

deltaPhi = 180/M;                      % Angular sampling rate phi
anglesPhi = 0:deltaPhi:180-deltaPhi;

gridSpacing =  -N/2:N/2;

%% Fully direct computations according to the adjoint formula
AdjointDFT_Point = 0;
for k = 1:K
        angleTheta = anglesTheta(k);
    for m = 1:M
         anglePhi = anglesPhi (m);
        for n = 1:N+1
            u = gridSpacing(n)*cosd(angleTheta)*cosd(anglePhi);
            v = gridSpacing(n)*cosd(angleTheta)*sind(anglePhi);
            w = gridSpacing(n)*sind(angleTheta);
%             if (1 == CheckPairsXZ(k,m,K,M))   %% This condition is for checking only XZ block
%           if (1 == CheckPairsXX(k,m,K,M))   %% This condition is for checking only XX block
%           if (angleTheta == 90)  %% This condition is for checking only for special line Z 
                AdjointDFT_Point = AdjointDFT_Point + SphericalGrid(k,m,n)*exp(+1i*2*pi*( r*u + c*v + d*w )/(N+1));        %% Observe the sign change as discussed in the paper for adjoint 
%             end
        end
    end
end

% %% Fully direct computations
% for l = 1:N+1            % x
%     for m = 1: N+1       % y
%         depthData(:) = I (l,m,:);
%         lineDataInner(m) = sum(depthData .* Map_n);
%     end
%     lineDataOuter(l) =  sum(lineDataInner .* Map_theta);
% end
% IDFT_Point  = sum( lineDataOuter .* Map_phi);
     
end

