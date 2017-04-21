clear all;
clc;
close all;


N = 8;   % Has to be even 
freqScale = 1; 2*pi/(N+1);
gridSpacing =  [-N/2:N/2]*freqScale;
gridMarkerSize = 6;
highlightMarkerSize = 8;

delay = .05;


%% Normal Grid
figure,axis equal
hold on
axis equal
axis ([-N/2 N/2 -N/2 N/2])
for j  = 1: N+1                             % Marking X-Y grid lines
    h =plot(gridSpacing(:),gridSpacing(j)*ones(1,N+1),'b');
    h =plot(gridSpacing(j)*ones(1,N+1),gridSpacing(:),'b');
end
for j  = 1: N+1                             % Plotting grid points
    pause(delay);
    h =plot(gridSpacing(j),gridSpacing(:),'bs');
    set(h,'MarkerFaceColor','g');
    set(h,'MarkerSize',gridMarkerSize)
end
h = plot(0,0,'bs');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize)

title('\bf{Square grid & Polar grid}') 
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
hold off

%% Plotting Polar grid
M = N + 2;                 % Total number of polar lines 
deltaTheta = 180/M;     % Angular sampling rate
angles = [0:deltaTheta:180];  
 
hold on 
for angle = angles 
    pause(delay) 
    h1= plot(gridSpacing*cosd(angle),gridSpacing*sind(angle),'ko');
    set(h1,'MarkerFaceColor','r');
    set(h1,'MarkerSize',gridMarkerSize) 
end
hold off  

 

%% First Scale sind(18)
alpha11 = cosd(deltaTheta);
figure,axis equal 
hold on 
axis equal
axis ([-N/2 N/2 -N/2 N/2])
for j  = 1: N+1                             % Marking X-Y grid lines
    h =plot(gridSpacing(:)*alpha11,gridSpacing(j)*ones(1,N+1),'b');
    h =plot(gridSpacing(j)*ones(1,N+1)*alpha11,gridSpacing(:),'b');
end
for j  = 1: N+1                             % Plotting grid points
    pause(delay);
    h =plot(gridSpacing(j)*alpha11,gridSpacing(:),'bs');
    set(h,'MarkerFaceColor','g');
    set(h,'MarkerSize',gridMarkerSize)
end
h = plot(0,0,'bs');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize)

title('\bf{1^{st} level scaling, x-axis by cos(\Delta\theta)}')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
hold off


hold on
for angle = angles
    pause(delay)
    h1= plot(gridSpacing*cosd(angle),gridSpacing*sind(angle),'ko');
    set(h1,'MarkerFaceColor','r');
    set(h1,'MarkerSize',gridMarkerSize)
end
hold off 

%% Next Scale in y direction
beta11 = sind(deltaTheta);
newspacing = abs(gridSpacing) * beta11;

newspacing(N/2+1) = 1;          % Scaling reset at zero

figure,axis equal
hold on 
axis equal
% axis ([-N/2 N/2 -N/2 N/2])
for j  = 1: N+1                             % Marking X-Y grid lines
    h =plot(gridSpacing(:)*alpha11,gridSpacing(j)*ones(1,N+1)*newspacing(j),'b');
    h =plot(gridSpacing(j)*ones(1,N+1)*alpha11,gridSpacing(:)*newspacing(j),'b');
end
for j  = 1: N+1                             % Plotting grid points
    pause(delay+.02);
    h =plot(gridSpacing(j)*alpha11,gridSpacing(:)*newspacing(j),'bs');
    set(h,'MarkerFaceColor','g');
    set(h,'MarkerSize',gridMarkerSize)
end
h = plot(0,0,'bs');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize)

title('\bf{1^{st} level scaling, x-axis by cos(\Delta\theta) then y axis by sin(\Delta\theta)}')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
hold off


hold on
for angle = angles
    pause(delay)
    h1= plot(gridSpacing*cosd(angle),gridSpacing*sind(angle),'ko');
    set(h1,'MarkerFaceColor','r');
    set(h1,'MarkerSize',gridMarkerSize)
end
hold off 
% axis tight



mywaitBar =waitbar(0,'Computing special grids');

%% Simple Rectangular Special Grid X-axis oriented
% beta11 = sind(deltaTheta);
newspacing = abs(gridSpacing) * beta11;
newspacing(N/2+1) = 1;          % Scaling reset at zero

figure,
% axis equal
hold on 
axis ([-6 6 -6 6])
for j  = 1: N+1                             % Marking X-Y grid lines
    if (gridSpacing(j) ~= 0)
    h =plot(gridSpacing(:)*alpha11,gridSpacing(j)*ones(1,N+1)*newspacing(j),'--b');
    h =plot(gridSpacing(j)*ones(1,N+1)*alpha11,gridSpacing(:)*newspacing(j),'--b');
    end
end
for j  = 1: N+1                             % Plotting grid points
    pause(delay+.02);
    if (gridSpacing(j) ~= 0)
    h =plot(gridSpacing(j)*alpha11,gridSpacing(:)*newspacing(j),'bs');
    set(h,'MarkerFaceColor','g');
    set(h,'MarkerSize',gridMarkerSize)
    end
end
h = plot(0,0,'bs');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize+2)

% title('\bf{The Butterfly X-Grid}')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
hold off

%% Special angles and radial scales
specialAngles1 = [ atan((gridSpacing (1:N/2).* newspacing(1:N/2)) ./ (gridSpacing(1:N/2)*alpha11))*180/pi ];
specialAngles1 = sort(specialAngles1);
% specialAngles2 = [ atan((gridSpacing (1:N/2).* newspacing(1:N/2)) ./ (gridSpacing(1:N/2)*cosd(180-deltaTheta)))*180/pi ];
specialAngles = [0 specialAngles1 90 -specialAngles1];
 
% Computing special scales
SpecialScales = [alpha11    alpha11./ cosd(abs(specialAngles1)) 0  alpha11./ cosd(abs(specialAngles1))]; 

 
hold on
for  i =1: length(specialAngles) 
    angle = specialAngles(i);
    scale = SpecialScales (i);
    pause(delay)
    h1 = plot(gridSpacing*cosd(angle)*scale,gridSpacing*sind(angle)*scale,'--b'); 
    h1 = plot(gridSpacing*cosd(angle)*scale,gridSpacing*sind(angle)*scale,'ko');
    set(h1,'MarkerFaceColor','r');
    set(h1,'MarkerSize',gridMarkerSize)
end
hold off 
 axis equal 
 
waitbar(.5);
 %% Simple Rectangular Special Grid X-axis oriented 
% beta11 = sind(deltaTheta);
newspacing = abs(gridSpacing) * beta11;
newspacing(N/2+1) = 1;          % Scaling reset at zero

figure,
% axis equal
hold on 
axis ([-6 6 -6 6])
for j  = 1: N+1                             % Marking X-Y grid lines
%     h =plot(gridSpacing(j)*ones(1,N+1)*newspacing(j), gridSpacing(:)*alpha11,'--b');
%     h =plot(gridSpacing(:)*newspacing(j),gridSpacing(j)*ones(1,N+1)*alpha11,'--b');
    h =plot(gridSpacing(j)*ones(1,N+1), gridSpacing(:),'--b');
    h =plot(gridSpacing(:)*newspacing(j),gridSpacing(j)*ones(1,N+1),'--b');
end
for j  = 1: N+1                             % Plotting grid points
    pause(delay+.02);
%     h =plot(gridSpacing(:)*newspacing(j),gridSpacing(j)*alpha11,'bs');
    h =plot(gridSpacing(:),gridSpacing(j),'bs');
    set(h,'MarkerFaceColor','g');
    set(h,'MarkerSize',gridMarkerSize)
end
h = plot(0,0,'bs');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize+2)

title('\bf{The Butterfly Y-Grid}')
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
hold off

%% Special angles and radial scales
specialAngles1 = [ atan((gridSpacing (1:N/2).* abs(gridSpacing(1:N/2)/ (N/2))) ./ (gridSpacing(1:N/2)*1))*180/pi ];
% specialAngles2 = [ atan((gridSpacing (1:N/2).* newspacing(1:N/2)) ./ (gridSpacing(1:N/2)*cosd(180-deltaTheta)))*180/pi ];
specialAngles = [0 specialAngles1 90 -specialAngles1];
 
% Computing special scales
SpecialScales = [1    1./ cosd(abs(specialAngles1)) 1   1./ cosd(abs(specialAngles1))]; 


hold on
for  i =1: length(specialAngles) 
    angle = specialAngles(i);
    scale = SpecialScales (i);
    pause(delay)
    h1 = plot(gridSpacing*sind(angle)*scale,gridSpacing*cosd(angle)*scale,'--b'); 
    h1 = plot(gridSpacing*sind(angle)*scale,gridSpacing*cosd(angle)*scale,'ko');
     set(h1,'MarkerFaceColor','r');
    set(h1,'MarkerSize',gridMarkerSize)
    
    h1 = plot(gridSpacing*cosd(angle)*scale,gridSpacing*sind(angle)*scale,'--b'); 
    h1 = plot(gridSpacing*cosd(angle)*scale,gridSpacing*sind(angle)*scale,'ko');
    
    set(h1,'MarkerFaceColor','r');
    set(h1,'MarkerSize',gridMarkerSize)
end
hold off 
axis tight
waitbar(1);
close(mywaitBar);