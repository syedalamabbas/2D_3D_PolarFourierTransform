clear all;
clc;
close all;


N = 8;   % Has to be even 
gridSpacing =  [-N/2:N/2];

% Figure specifications
gridMarkerSize = 3;
highlightMarkerSize = 8;
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
M = 10;                 % Total number of polar lines
deltaTheta = 180/M;     % Angular sampling rate
deltaPhi  = deltaTheta;

gridSpacing = -4:4;
angles = [0:deltaTheta:180-deltaTheta];  
anglesPhi = [0:deltaPhi:180-deltaTheta];
count = 0;
figure,
rotate3d on;
hold on
for anglePhi = anglesPhi
    for angleTheta = angles
%           pause(delay+30* .05)
        h1 =plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'g', 'LineWidth', 1);
        
        if (count == 0)
            rotate(h1,[1 1 1],[ -13]);
        end 
%             pause(delay+3* .05)
        h1= plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'ko');
        
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


%% Plotting single slice X-Y , perpendicular to Z axis
figure,
hold on
for anglePhi = anglesPhi 
    for angleTheta = angles
        pause(delay)
        if (anglePhi ==  18)
            h1 =plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'g', 'LineWidth', 1);
            h1= plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'ko');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerSize',gridMarkerSize+1)
        end
    end
end
title('\bf{Single slice of a spherical grid at \phi = 18^0 }')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
zlabel('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on
hold off
%% Plotting Spherical grid  X axis oriented block at 18 degrees apart 

figure,
hold on
for anglePhi = anglesPhi 
    for angleTheta = angles 
        pause(delay)
        if (anglePhi == 18 || anglePhi == (180-18)  )
            h1 =plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'g', 'LineWidth', 1);
            h1= plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'ko');
%             h1 =plot3(gridSpacing*sind(anglePhi)*cosd(angleTheta),gridSpacing*sind(anglePhi)*sind(angleTheta),gridSpacing*cosd(anglePhi),'g', 'LineWidth', 1);
%             h1= plot3(gridSpacing*sind(anglePhi)*cosd(angleTheta),gridSpacing*sind(anglePhi)*sind(angleTheta),gridSpacing*cosd(anglePhi),'ko');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerSize',gridMarkerSize+1)
        end 
    end 
end 

Px = 3.618;
Py = 1.176;
Pz = 1.236;

P1 = [ Px , Py, Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, -Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , -Py, Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, Pz];
P2 = [Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [-Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [-Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , -Py, Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , -Py, -Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

title('\bf{Spherical grid X axis oriented block at angle \theta = 18^0 , \phi = [18^0],[162^0]}')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}') 
zlabel('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on
hold off  




%% Plotting Spherical grid  X axis oriented block at 18 degrees apart but at 2nd angle

sliceAtangle = 18; 
figure,
hold on
for anglePhi = anglesPhi
    for angleTheta = angles 
        pause(delay)
        if (anglePhi == sliceAtangle || anglePhi == (180-sliceAtangle)  )
            h1 =plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'g', 'LineWidth', 1);
            h1= plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'ko');
%                h1 =plot3(gridSpacing*sind(anglePhi)*cosd(angleTheta),gridSpacing*sind(anglePhi)*sind(angleTheta),gridSpacing*cosd(anglePhi),'g', 'LineWidth', 1);
%             h1= plot3(gridSpacing*sind(anglePhi)*cosd(angleTheta),gridSpacing*sind(anglePhi)*sind(angleTheta),gridSpacing*cosd(anglePhi),'ko');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerSize',gridMarkerSize+1)
        end
    end
end

Px = 3.078;
Py = 1;
Pz = 2.351;

P1 = [ Px , Py, Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, -Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , -Py, Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, Pz];
P2 = [Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [-Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [-Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , -Py, Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , -Py, -Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

title('\bf{Spherical grid X axis oriented block, 2^{nd} angle, \theta = 36^0 , \phi = [18^0],[162^0] }')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}') 
zlabel('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on
hold off  


%% Plotting Spherical grid  Level 1 Z axis oriented block at 18 degrees apart at 1st angle

sliceAtangle = 18;
figure,
hold on
for anglePhi = anglesPhi
    for angleTheta = angles 
        pause(delay)
        if (anglePhi == sliceAtangle || anglePhi == (180-sliceAtangle)  )
            h1 =plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'g', 'LineWidth', 1);
            h1= plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'ko');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerSize',gridMarkerSize+1)
        end
    end
end

Px = 1.176;
Py = .382;
Pz = 3.804;

P1 = [ Px , -Py, -Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, -Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, Pz];
P2 = [-Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [-Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , -Py, Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , -Py, -Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

title('\bf{Spherical grid Z axis oriented block, 1^{st} angle, \theta = 18^0, \phi = (18^0,162^0) }')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}') 
zlabel('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on
hold off  


%% Plotting Spherical grid  Level 1 Z axis oriented block at 18 degrees apart at 2nd angle

sliceAtangle = 18; 
figure,
hold on
for anglePhi = anglesPhi
    for angleTheta = angles 
        pause(delay)
        if (anglePhi == sliceAtangle || anglePhi == (180-sliceAtangle)  )
            h1 =plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'g', 'LineWidth', 1);
            h1= plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'ko');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerSize',gridMarkerSize+1)
        end
    end
end

Px = 2.236;
Py = .7265;
Pz = 3.236;

P1 = [ Px , -Py, -Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, -Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, Pz];
P2 = [-Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [-Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , -Py, Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , -Py, -Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

title('\bf{Spherical grid Z axis oriented block, 2^{nd} angle, \theta = 36^0 , \phi = [18^0],[162^0] }')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}') 
zlabel('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on
hold off  


%% Plotting Spherical grid   Y axis oriented block
figure,
hold on
for anglePhi = anglesPhi
    for angleTheta = angles
        pause(delay)
        if (anglePhi == (90 - 18) || anglePhi == (90+18)  )
            h1 =plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'g', 'LineWidth', 1);
            h1= plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'ko');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerSize',gridMarkerSize+1)
        end
    end
end


Px = 1.176 ;
Py = 3.618;
Pz = 1.236;

P1 = [ Px , Py, Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ Px , -Py, -Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , -Py, Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, Pz];
P2 = [-Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ Px , -Py, Pz];
P2 = [Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , Py, Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , Py, Pz];
P2 = [-Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , -Py, Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , -Py, -Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , Py, Pz];
P2 = [Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


title('\bf{Spherical grid  Y axis oriented block at angle \theta = 18^0}')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}') 
zlabel('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on
hold off  


%% Plotting Spherical grid   Y axis oriented block, 2nd angle
figure,
hold on
for anglePhi = anglesPhi
    for angleTheta = angles
        pause(delay)
        if (anglePhi == (90 - 18) || anglePhi == (90+18)  )
            h1 =plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'g', 'LineWidth', 1);
            h1= plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'ko');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerSize',gridMarkerSize+1)
        end
    end
end


Px = 1 ;
Py = 3.078;
Pz = 2.351;

P1 = [ Px , Py, Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ Px , -Py, -Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , -Py, Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, Pz];
P2 = [-Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ Px , -Py, Pz];
P2 = [Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , Py, Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , Py, Pz];
P2 = [-Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , -Py, Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , -Py, -Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , Py, Pz];
P2 = [Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


title('\bf{Spherical grid  Y axis oriented block, angle \theta = 36^0}')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}') 
zlabel('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on
hold off  


%% Plotting Spherical grid  Level 1 YZ axis oriented block at 90 +- 18 degrees apart at 1st angle

sliceAtangle = 18;
figure,
hold on
for anglePhi = anglesPhi
    for angleTheta = angles 
        pause(delay)
        if (anglePhi == 90-sliceAtangle || anglePhi == (90+sliceAtangle)  )
            h1 =plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'g', 'LineWidth', 1);
            h1= plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'ko');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerSize',gridMarkerSize+1)
        end
    end
end

Px = .382;
Py = 1.176;
Pz = 3.804;

P1 = [ Px , -Py, -Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, -Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, Pz];
P2 = [-Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [-Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , -Py, Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , -Py, -Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

title('\bf{Spherical grid Y-Z axis oriented block, 1^{st} angle, \theta = 90 +- 18 , \phi = [18^0],[162^0] }')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}') 
zlabel('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on
hold off  


%% Plotting Spherical grid  Level 1 Z axis oriented block at 18 degrees apart at 2nd angle

sliceAtangle = 2* 18;
figure,
hold on
for anglePhi = anglesPhi
    for angleTheta = angles 
        pause(delay)
        if (anglePhi == 90 - sliceAtangle || anglePhi == (90+sliceAtangle)  )
            h1 =plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'g', 'LineWidth', 1);
            h1= plot3(gridSpacing*sind(angleTheta)*cosd(anglePhi),gridSpacing*sind(angleTheta)*sind(anglePhi),gridSpacing*cosd(angleTheta),'ko');
            set(h1,'MarkerFaceColor','r');
            set(h1,'MarkerSize',gridMarkerSize+1)
        end
    end
end

Px = 1.382;
Py = 1.902;
Pz = 3.236;

P1 = [ Px , -Py, -Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, -Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ Px , -Py, Pz];
P2 = [-Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [-Px, -Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)


P1 = [ -Px , -Py, Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , -Py, -Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [-Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [Px, Py, Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, Pz];
P2 = [Px, -Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

P1 = [ -Px , Py, -Pz];
P2 = [Px, Py, -Pz];
pts = [P1; P2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'LineWidth', 1.7)

title('\bf{Spherical grid Y-Z axis oriented block, 2^{nd} angle, \theta = 90 +- 36 , \phi = [18^0],[162^0] }')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}') 
zlabel('\bf{z \rightarrow}')
axis equal
axis ([-N/2 N/2 -N/2 N/2 -N/2 N/2])
grid on
hold off  

