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
% for j  = 1: N+1                             % Plotting grid points
%     pause(delay);
%     h =plot(gridSpacing(j),gridSpacing(:),'bs');
%     set(h,'MarkerFaceColor','g');
%     set(h,'MarkerSize',gridMarkerSize)
% end
% h = plot(0,0,'bs');                 % DC value
% set(h,'MarkerFaceColor','k');
% set(h,'MarkerSize',gridMarkerSize) 

% title('\bf{Square grid & Polar grid & Spiral grid}') 
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
hold off
 


%% Plotting Polar grid
M = 10;                 % Total number of polar lines 
deltaTheta = 180/M;     % Angular sampling rate
angles = [0:deltaTheta:180];  
 
hold on 
for angle = angles 
    pause(delay) 
    h1= plot(gridSpacing*cosd(angle),gridSpacing*sind(angle),'ko');
    set(h1,'MarkerFaceColor','g');
    set(h1,'MarkerSize',gridMarkerSize+1) 
end
hold off  



%% Plotting Spiral grid
M = 10;                 % Total number of polar lines 
deltaTheta = pi/M;     % Angular sampling rate
angles = [0:deltaTheta:10*pi];  
 
offset = (0:1:11)*pi;
count = 1;
k = .13;   % constant
hold on 
for angle = angles 
    pause(delay)  
    
    h1= plot(k *angle*cos(angle),k *angle*sin(angle),'ko');
    set(h1,'MarkerFaceColor','r');
    set(h1,'MarkerSize',gridMarkerSize) 
end

Newspacing = deltaTheta + offset;
plot(k *Newspacing.*cos(Newspacing),k *Newspacing.*sin(Newspacing),'--c');
  
Newspacing = 9*deltaTheta + offset;
plot(k *Newspacing.*cos(Newspacing),k *Newspacing.*sin(Newspacing),'--c');
hold off  