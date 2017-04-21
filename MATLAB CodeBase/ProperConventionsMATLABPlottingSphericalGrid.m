clear all;
clc;
close all;


N = 8;   % Has to be even
gridSpacing =  [-N/2:N/2];

% Figure specifications
gridMarkerSize = 3;
% highlightMarkerSize = 8;
delay = .03;


%% Normal Grid
figure,axis equal
hold on
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
for k =1 :N+1
    for j  = 1: N+1                             % Marking X-Y grid lines
        h =plot3(gridSpacing(k)*ones(1,N+1),gridSpacing(:),gridSpacing(j)*ones(1,N+1),'b', 'LineWidth', .2);
        h =plot3(gridSpacing(k)*ones(1,N+1),gridSpacing(j)*ones(1,N+1),gridSpacing(:),'b', 'LineWidth', .2);
        h =plot3(gridSpacing(:),gridSpacing(k)*ones(1,N+1),gridSpacing(j)*ones(1,N+1),'b', 'LineWidth', .2);
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
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
%% Plotting transparent spheres
shadeVal = .5;
points = 30;
colorVec = [1 1 1];

for l =0; 0:3
r = N/2-l;
[x,y,z] = sphere(points);
x = x*r ;
y = y*r ;
z = z*r ;
lightGrey = shadeVal*colorVec; % It looks better if the lines are lighter
surface(x,y,z,'FaceColor', 'none','EdgeColor',lightGrey)

end


for anglePhi = anglesPhi  % As defined in the slides , phi (Azimuth or latitude ?) first then  theta (elevation or zenith) for spherical divisions
    for angleTheta = angles
%          pause(delay+30* .05)
        h1 =plot3(gridSpacing*cosd(angleTheta)*cosd(anglePhi),gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta),'g', 'LineWidth', 1);  % as defined in the paper
        
%                      pause(delay+3* .05)
        h1= plot3(gridSpacing*cosd(angleTheta)*cosd(anglePhi),gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta),'ko');
        
        
        set(h1,'MarkerFaceColor','r');
        set(h1,'MarkerSize',gridMarkerSize+2)
        count = count +1;
    end
end

title('\bf{Spherical grid}')
% title('\bf{Spherical grid top view}')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
zlabel('\bf{z \rightarrow}')

grid on
hold off


%% Just a sphere
r = N/2;

[x,y,z] = sphere(points);
x = x*r ;
y = y*r ;
z = z*r ;
lightGrey = shadeVal*colorVec; % It looks better if the lines are lighter
figure, surface(x,y,z,'FaceColor', 'none','EdgeColor',lightGrey)
title('\bf{Spherical grid}')
% title('\bf{Spherical grid top view}')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
zlabel('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on

%% Plotting single slice X-Y at an angle phi = 18, perpendicular to Z axis
currentAnglePhi = deltaPhi;
figure,
title(['\bf{Single polar slice of a spherical grid at \phi = }', num2str(currentAnglePhi), '^0'])
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
zlabel ('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on
hold on
for anglePhi = anglesPhi
    for angleTheta = angles
        
        pause(delay)
        if ( anglePhi == currentAnglePhi) % angleTheta == currentAngleTheta ) %anglePhi
            %             h1 =plot3(gridSpacing*cosd(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(anglePhi),'g', 'LineWidth', 1);
            %             h1= plot3(gridSpacing*cosd(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(anglePhi),'ko');
            
            h1 =plot3(gridSpacing*cosd(angleTheta)*cosd(anglePhi),gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta),'g', 'LineWidth', 1);
            h1= plot3(gridSpacing*cosd(angleTheta)*cosd(anglePhi),gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta),'ko');
            
%             h1 =plot3(gridSpacing*cosd(anglePhi)*cosd(angleTheta),gridSpacing*sind(anglePhi)*cosd(angleTheta),gridSpacing*sind(angleTheta),'g', 'LineWidth', 1);
%             h1= plot3(gridSpacing*cosd(anglePhi)*cosd(angleTheta),gridSpacing*sind(anglePhi)*cosd(angleTheta),gridSpacing*sind(angleTheta),'ko');
            
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerSize',gridMarkerSize+1)
        end
    end
end

hold off

%% Designing the Angles
angleTheta_multiple = 2;   % Elevation or Zenith angles 
anglePhi_multiple = 2;     % Latitude or Azimuth angles  

currentAngleThetaWithFactor = angleTheta_multiple*deltaTheta; %angleTheta can only be changed
currentAnglePhiWithFactor = anglePhi_multiple* deltaPhi;


%% Plotting Spherical grid  X axis oriented block at 18 degrees apart
% This is set of YZ blocks tiled along X -axis

phi1 =  currentAnglePhiWithFactor;         % This is how XX block looks like
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
for anglePhi = anglesPhi
    for angleTheta = angles
        pause(delay)
        if (anglePhi == phi1 || anglePhi == phi2  )
            h1 =plot3(gridSpacing*cosd(anglePhi)*cosd(angleTheta),gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta),'g', 'LineWidth', 1);
            h1= plot3(gridSpacing*cosd(anglePhi)*cosd(angleTheta),gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta),'ko');
            
%             h1 =plot3(gridSpacing*cosd(anglePhi)*cosd(angleTheta),gridSpacing*cosd(anglePhi)*sind(angleTheta),gridSpacing*sind(anglePhi),'g', 'LineWidth', 1);
%             h1= plot3(gridSpacing*cosd(anglePhi)*cosd(angleTheta),gridSpacing*cosd(anglePhi)*sind(angleTheta),gridSpacing*sind(anglePhi),'ko');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerSize',gridMarkerSize+1)
        end
    end
end

alpha_factor = gridSpacing*cosd(currentAnglePhiWithFactor)*cosd(currentAngleThetaWithFactor)  ;       % Scaling needed in X-axis as defined in the schematic of the paper
beta_factor  = gridSpacing*sind(currentAnglePhiWithFactor)*cosd(currentAngleThetaWithFactor);
gamma_factor = gridSpacing*sind(currentAngleThetaWithFactor);

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
 

%% Plotting Spherical grid  X axis oriented block at 18 degrees apart
% This is set of YZ blocks tiled along X -axis

% currentAngleThetaWithFactor = 18; %angleTheta can only be changed
% currentAnglePhiWithFactor = angleTheta_multiple*18;
phi1 =  currentAnglePhiWithFactor;
phi2 = 180 - currentAnglePhiWithFactor;
theta1 = 90 - currentAngleThetaWithFactor;
theta2 = 90 + currentAngleThetaWithFactor;

figure,
title(['\bf{Spherical grid XZ axis oriented block at angle \theta = (', num2str(theta1),'^0,', num2str(theta2),'^0), \phi =(', num2str(phi1), '^0,', num2str(phi2), '^0)}'])
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
zlabel('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on
hold on
for anglePhi = anglesPhi
    for angleTheta = angles
        pause(delay)
        if (anglePhi == phi1 || anglePhi == phi2 )
            h1 =plot3(gridSpacing*cosd(anglePhi)*cosd(angleTheta),gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta),'g', 'LineWidth', 1);
            h1= plot3(gridSpacing*cosd(anglePhi)*cosd(angleTheta),gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta),'ko');
            
%             h1 =plot3(gridSpacing*cosd(anglePhi)*cosd(angleTheta),gridSpacing*cosd(anglePhi)*sind(angleTheta),gridSpacing*sind(anglePhi),'g', 'LineWidth', 1);
%             h1= plot3(gridSpacing*cosd(anglePhi)*cosd(angleTheta),gridSpacing*cosd(anglePhi)*sind(angleTheta),gridSpacing*sind(anglePhi),'ko');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerSize',gridMarkerSize+1)
        end
    end
end

alpha_factor = gridSpacing*sind(currentAngleThetaWithFactor)*cosd(currentAnglePhiWithFactor);        % Scaling needed in X-axis
beta_factor  = gridSpacing*sind(currentAngleThetaWithFactor)*sind(currentAnglePhiWithFactor);
gamma_factor = gridSpacing*cosd(currentAngleThetaWithFactor)  ;

plot3(alpha_factor, beta_factor, gamma_factor, 'LineWidth', 1.7)
plot3(-alpha_factor, beta_factor, -gamma_factor, 'LineWidth', 1.7)

for i=1: N+1
    P1 = [alpha_factor(1, i) , beta_factor(1,i) , gamma_factor(1,i)];
    P2 = [alpha_factor(1,i),  -beta_factor(1,i) , gamma_factor(1,i) ];
    pts = [P1; P2];
    plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)
end
for i=1: N+1
    P1 = [-alpha_factor(1, i) , beta_factor(1,i) , gamma_factor(1,i)];
    P2 = [alpha_factor(1,i), beta_factor(1,i) , gamma_factor(1,i) ];
    pts = [P1; P2];
    plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)
end

plot3(alpha_factor, beta_factor, -gamma_factor, 'LineWidth', 1.7)
plot3(alpha_factor, -beta_factor, -gamma_factor, 'LineWidth', 1.7)

for i=1: N+1
    P1 = [alpha_factor(1, i) , beta_factor(1,i) , -gamma_factor(1,i)];
    P2 = [alpha_factor(1,i),  -beta_factor(1,i) , -gamma_factor(1,i) ];
    pts = [P1; P2];
    plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)
end
for i=1: N+1
    P1 = [-alpha_factor(1, i) , beta_factor(1,i) , -gamma_factor(1,i)];
    P2 = [alpha_factor(1,i), beta_factor(1,i) , -gamma_factor(1,i) ];
    pts = [P1; P2];
    plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)
end
hold off

%% Plotting Spherical grid  X axis oriented block at 18 degrees apart
% This is set of YZ blocks tiled along X -axis

% currentAngleThetaWithFactor = angleTheta_multiple*18; %angleTheta can only be changed
% currentAnglePhiWithFactor =  18;

phi1 = 90 -  currentAnglePhiWithFactor;  
phi2 = 90 +  currentAnglePhiWithFactor;
theta1 = currentAngleThetaWithFactor;
theta2 = 180 - currentAngleThetaWithFactor;

figure,
title(['\bf{Spherical grid YY axis oriented block at angle \theta = (', num2str(theta1),'^0,', num2str(theta2),'^0), \phi =(', num2str(phi1), '^0,', num2str(phi2), '^0)}'])
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
zlabel('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on
hold on
for anglePhi = anglesPhi
    for angleTheta = angles
        pause(delay)
        if (anglePhi == phi1 || floor( anglePhi) == floor(phi2) )
            h1 =plot3(gridSpacing*cosd(angleTheta)*cosd(anglePhi),gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta),'g', 'LineWidth', 1);
            h1= plot3(gridSpacing*cosd(angleTheta)*cosd(anglePhi),gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta),'ko');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerSize',gridMarkerSize+1)
        end
    end
end


alpha_factor = gridSpacing*sind(currentAnglePhiWithFactor)*cosd(currentAngleThetaWithFactor);% Scaling needed in X-axis
beta_factor  = gridSpacing*cosd(currentAnglePhiWithFactor)*cosd(currentAngleThetaWithFactor);
gamma_factor = gridSpacing*sind(currentAngleThetaWithFactor);

plot3(alpha_factor, beta_factor, gamma_factor, 'LineWidth', 1.7)
plot3(-alpha_factor, beta_factor, -gamma_factor, 'LineWidth', 1.7)

for i=1: N+1
    P1 = [-alpha_factor(1, i) , beta_factor(1,i) , gamma_factor(1,i)];
    P2 = [alpha_factor(1,i),  beta_factor(1,i) , gamma_factor(1,i) ];
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
    P1 = [alpha_factor(1, i) ,- beta_factor(1,i) , gamma_factor(1,i)];
    P2 = [alpha_factor(1,i),  -beta_factor(1,i) , -gamma_factor(1,i) ];
    pts = [P1; P2];
    plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)
end
for i=1: N+1
    P1 = [alpha_factor(1, i) ,- beta_factor(1,i) , gamma_factor(1,i)];
    P2 = [-alpha_factor(1,i), - beta_factor(1,i) , gamma_factor(1,i) ];
    pts = [P1; P2];
    plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)
end
hold off


%% Plotting Spherical grid  X axis oriented block at 18 degrees apart
% This is set of YZ blocks tiled along X -axis

% currentAngleThetaWithFactor = 18; %angleTheta can only be changed
% currentAnglePhiWithFactor =  angleTheta_multiple*18;
phi1 = 90 - currentAnglePhiWithFactor;
phi2 = 90 + currentAnglePhiWithFactor;
theta1 = 90 - currentAngleThetaWithFactor;
theta2 = 90 + currentAngleThetaWithFactor;

figure,
title(['\bf{Spherical grid YZ axis oriented block at angle \theta = (', num2str(theta1),'^0,', num2str(theta2),'^0), \phi =(', num2str(phi1), '^0,', num2str(phi2), '^0)}'])
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
zlabel('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on
hold on
for anglePhi = anglesPhi
    for angleTheta = angles
        pause(delay)
        if (anglePhi == phi1 || floor(anglePhi) == floor(phi2) )
            h1 =plot3(gridSpacing*cosd(angleTheta)*cosd(anglePhi),gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta),'g', 'LineWidth', 1);
            h1= plot3(gridSpacing*cosd(angleTheta)*cosd(anglePhi),gridSpacing*cosd(angleTheta)*sind(anglePhi),gridSpacing*sind(angleTheta),'ko');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerSize',gridMarkerSize+1)
        end
    end
end


alpha_factor = gridSpacing*sind(currentAngleThetaWithFactor)*sind(currentAnglePhiWithFactor);      % Scaling needed in X-axis
beta_factor  = gridSpacing*sind(currentAngleThetaWithFactor)*cosd(currentAnglePhiWithFactor);
gamma_factor = gridSpacing*cosd(currentAngleThetaWithFactor)  ;

plot3(alpha_factor, beta_factor, gamma_factor, 'LineWidth', 1.7)
plot3(-alpha_factor, beta_factor, -gamma_factor, 'LineWidth', 1.7)

for i=1: N+1
    P1 = [alpha_factor(1, i) , beta_factor(1,i) , gamma_factor(1,i)];
    P2 = [alpha_factor(1,i),  -beta_factor(1,i) , gamma_factor(1,i) ];
    pts = [P1; P2];
    plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)
end
for i=1: N+1
    P1 = [-alpha_factor(1, i) , beta_factor(1,i) , gamma_factor(1,i)];
    P2 = [alpha_factor(1,i), beta_factor(1,i) , gamma_factor(1,i) ];
    pts = [P1; P2];
    plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)
end

plot3(alpha_factor, beta_factor, -gamma_factor, 'LineWidth', 1.7)
plot3(alpha_factor, -beta_factor, -gamma_factor, 'LineWidth', 1.7)

for i=1: N+1
    P1 = [alpha_factor(1, i) , beta_factor(1,i) , -gamma_factor(1,i)];
    P2 = [alpha_factor(1,i),  -beta_factor(1,i) , -gamma_factor(1,i) ];
    pts = [P1; P2];
    plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)
end
for i=1: N+1
    P1 = [-alpha_factor(1, i) , beta_factor(1,i) , -gamma_factor(1,i)];
    P2 = [alpha_factor(1,i), beta_factor(1,i) , -gamma_factor(1,i) ];
    pts = [P1; P2];
    plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)
end
hold off