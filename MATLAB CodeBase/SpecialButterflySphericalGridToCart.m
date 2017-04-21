clear all;
clc;
close all;


N = 8;   % Has to be even
gridSpacing =  [-N/2:N/2];

% Figure specifications
gridMarkerSize = 4;
% highlightMarkerSize = 8;
delay = .03;


%% Normal Grid
figure,axis equal
hold on
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
for k =1 :N+1 
    for j  = 1: N+1                             % Marking X-Y grid lines
        h =plot3(gridSpacing(k)*ones(1,N+1),gridSpacing(:),gridSpacing(j)*ones(1,N+1),'--b', 'LineWidth', .2);
        h =plot3(gridSpacing(k)*ones(1,N+1),gridSpacing(j)*ones(1,N+1),gridSpacing(:),'--b', 'LineWidth', .2);
        h =plot3(gridSpacing(:),gridSpacing(k)*ones(1,N+1),gridSpacing(j)*ones(1,N+1),'--b', 'LineWidth', .2);
    end
end
for k =1 :N+1
    for j  = 1: N+1                             % Plotting grid points
        %         pause(delay);
        h =plot3(gridSpacing(k),gridSpacing(j),gridSpacing(:),'bs');
        set(h,'MarkerFaceColor','g');
        set(h,'MarkerSize',gridMarkerSize)
    end
end
h = plot(0,0,'bs');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize)

title('\bf{Cube grid}')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
zlabel('\bf{z \rightarrow}')
hold off


%% Plotting Simple Spherical grid
M = N+2;                 % Total number of polar lines
deltaTheta = 180/M;     % Angular sampling rate
deltaPhi  = deltaTheta;

gridSpacing = -N/2:N/2;
angles = [0:deltaTheta:180-deltaTheta];
anglesPhi = [0:deltaPhi:180-deltaTheta];
count = 0;
figure,
rotate3d on;
hold on
for anglePhi = anglesPhi
    for angleTheta = angles
        %           pause(delay+30* .05)
        h1 =plot3(gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(anglePhi),'g', 'LineWidth', 1);
        
        if (count == 0)
            rotate(h1,[1 1 1],[ -13]);
        end
        %             pause(delay+3* .05)
        h1= plot3(gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(anglePhi),'ko');
        
        set(h1,'MarkerFaceColor','r');
        set(h1,'MarkerSize',gridMarkerSize+1)
        count = count +1;
    end
end

title('\bf{Spherical grid top view}')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
zlabel('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on
hold off

%% Designing the Angles
angleTheta_multiple = 1;
anglePhi_multiple = 1;

currentAngleThetaWithFactor = angleTheta_multiple*deltaTheta; %angleTheta can only be changed
currentAnglePhiWithFactor = anglePhi_multiple* deltaPhi;


%% Plotting Spherical grid  X axis oriented block at 18 degrees apart
% This is set of YZ blocks tiled along X -axis

phi1 =  currentAnglePhiWithFactor;
phi2 = 180 - currentAnglePhiWithFactor;
theta1 =  currentAngleThetaWithFactor;
theta2 = 180 - currentAngleThetaWithFactor;

figure,
title(['\bf{Spherical grid XX axis oriented block at angle \theta = (', num2str(theta1),'^0,', num2str(theta2),'^0), \phi =(', num2str(phi1), '^0,', num2str(phi2), '^0)}'])
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
zlabel('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on
hold on
for angleTheta = angles
    for anglePhi = anglesPhi
        pause(delay)
        if (angleTheta == theta1 || angleTheta == theta2  )
            h1 =plot3(gridSpacing*cosd(anglePhi)*cosd(angleTheta),gridSpacing*cosd(anglePhi)*sind(angleTheta),gridSpacing*sind(anglePhi),'g', 'LineWidth', 1);
            h1= plot3(gridSpacing*cosd(anglePhi)*cosd(angleTheta),gridSpacing*cosd(anglePhi)*sind(angleTheta),gridSpacing*sind(anglePhi),'ko');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerSize',gridMarkerSize+1)
        end
    end
end

alpha = cosd(currentAngleThetaWithFactor)*cosd(currentAnglePhiWithFactor) ;
beta = sind(currentAngleThetaWithFactor)*cosd(currentAnglePhiWithFactor);
gamma = sind(currentAnglePhiWithFactor);

alpha_factor = gridSpacing*alpha  ;       % Scaling needed in X-axis
beta_factor  = gridSpacing*beta ;
gamma_factor = gridSpacing*gamma ;

plot3(alpha_factor, beta_factor, gamma_factor, 'LineWidth', 1.7)
plot3(-alpha_factor, beta_factor, -gamma_factor, 'LineWidth', 1.7)

for i=1: N+1
    P1 = [alpha_factor(1, i) , beta_factor(1,i) , gamma_factor(1,i)];
    P2 = [alpha_factor(1,i),  -beta_factor(1,i) , gamma_factor(1,i) ];
    pts = [P1; P2];
    plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)
end
for i=1: N+1
    P1 = [alpha_factor(1, i) , beta_factor(1,i) , gamma_factor(1,i)];
    P2 = [alpha_factor(1,i), beta_factor(1,i) , -gamma_factor(1,i) ];
    pts = [P1; P2];
    plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)
end

plot3(alpha_factor, beta_factor, -gamma_factor, 'LineWidth', 1.7)
plot3(alpha_factor, -beta_factor, -gamma_factor, 'LineWidth', 1.7)

for i=1: N+1
    P1 = [-alpha_factor(1, i) , beta_factor(1,i) , gamma_factor(1,i)];
    P2 = [-alpha_factor(1,i),  beta_factor(1,i) , -gamma_factor(1,i) ];
    pts = [P1; P2];
    plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)
end
for i=1: N+1
    P1 = [-alpha_factor(1, i) , beta_factor(1,i) , gamma_factor(1,i)];
    P2 = [-alpha_factor(1,i),  -beta_factor(1,i) , gamma_factor(1,i) ];
    pts = [P1; P2];
    plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)
end
hold off

%% Butterfly Grid


figure,axis equal
h = plot(0,0,'bs');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize + 2)
grid on
% axis tight
% title('\bf{Butterfly grid}')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}') 
zlabel('\bf{z \rightarrow}')
hold on
axis equal
% axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
for k = 1 :N+1
    for j  = 1: N+1                             % Marking X-Y grid lines
        h =plot3(gridSpacing(k)*ones(1,N+1)*alpha,gridSpacing* abs(gridSpacing(k))*beta,gridSpacing(j)*ones(1,N+1)* abs(gridSpacing(k))*gamma,'--b', 'LineWidth', .2);
        h =plot3(gridSpacing(k)*ones(1,N+1)*alpha,gridSpacing(j)*ones(1,N+1)* abs(gridSpacing(k))*beta,gridSpacing* abs(gridSpacing(k))*gamma,'--b', 'LineWidth', .2);
        %         h =plot3(gridSpacing(:)*alpha,gridSpacing(k)*ones(1,N+1)* abs(gridSpacing(k))*beta,gridSpacing(j)*ones(1,N+1)* abs(gridSpacing(k))*gamma,'--b', 'LineWidth', .2);
    end
end
for k = 1:N+1 
    for j  = 1: N+1                             % Plotting grid points
        %         pause(delay);
        h =plot3(gridSpacing(k)*alpha,gridSpacing(j)* abs(gridSpacing(k))*beta,gridSpacing(:)* abs(gridSpacing(k))*gamma,'ks');
        set(h,'MarkerFaceColor','g');
        set(h,'MarkerSize',gridMarkerSize+1)
    end
end

%% Special computations of angles
specialAngles1 = [ atan((gridSpacing (1:N/2).* abs(gridSpacing(1:N/2))*beta ) ./ (gridSpacing(1:N/2)*alpha))*180/pi ];
specialAnglesTheta = [ specialAngles1 0 -fliplr(specialAngles1)];

x_side = gridSpacing(N/2)*ones(1,N+1)*alpha;
y_side = gridSpacing(1:N+1)* abs(gridSpacing(N/2))*beta;
z_side = (gridSpacing (1:N+1).* abs(gridSpacing(N/2))*gamma );
adjSide =  sqrt( x_side.^2 + y_side.^2);

specialAnglesPhi = zeros(N+1, N+1);
specialScales = zeros(N+1, N+1);

for k = 1: N+1
    specialAngles1 = [ atan( z_side ./   adjSide(k) )*180/pi  ];      % computing the angle subtended by line w.r.t X-Y plane
    specialAnglesPhi(k,: )  =  specialAngles1;
    specialScales(k,:) = sqrt( z_side.^2 + adjSide(k)^2);
end

 
for counter = 1: length(specialAnglesTheta )
    angleTheta = specialAnglesTheta(counter) ;
    specialAnglesPhiSingle = specialAnglesPhi(counter, :);
    
    SpecialScales = specialScales(counter, :);
        
    for k=1: length(specialAnglesPhiSingle)
        anglePhi = specialAnglesPhiSingle(k);
        scale =  SpecialScales(k);
        pause(delay)
        h1 =plot3(scale*gridSpacing*cosd(anglePhi)*cosd(angleTheta),scale*gridSpacing*cosd(anglePhi)*sind(angleTheta),scale*gridSpacing*sind(anglePhi),'-.r','LineWidth', .4);
    end
end
hold off
view(-5,46)
pause(3)

%  rotate3d on;
%  direction = [1 1 0];
%  rotate(h1,direction, 25)

% for fancy = 1: 1000
%     pause(delay)
%      camorbit(2,2,'camera')
%      drawnow
% end
