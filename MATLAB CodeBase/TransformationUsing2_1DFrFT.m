close all;
clear all;
clc;

%% First Spatial grid
N = 8;            % even number of data points used  
freqScale = 1; 2*pi/(N+1);
gridSpacing =  [-N/2:N/2]*freqScale;
gridMarkerSize = 6;
highlightMarkerSize = 8;

figure,
hold on
axis equal
axis ([-N/2 N/2 -N/2 N/2+.002])
for j  = 1: N+1                             % Marking X-Y grid lines
    h =plot(gridSpacing(:),gridSpacing(j)*ones(1,N+1),'b');
    h =plot(gridSpacing(j)*ones(1,N+1),gridSpacing(:),'b');
end
for j  = 1: N+1                             % Plotting grid points
    pause(.03);
    h =plot(gridSpacing(j),gridSpacing(:),'bo');
    set(h,'MarkerFaceColor','g');
    set(h,'MarkerSize',gridMarkerSize)
end
h = plot(0,0,'bo');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize)

title('Original Image:(N+1) \times (N+1), N = 4 is even')
xlabel ('x \rightarrow')
ylabel ('y \rightarrow')
hold off


%% Pick an arbitrary point on the grid
%   desiredPoint = [ 0,0];
desiredPoint = [ 1.24 , 2.56 ]*freqScale;
% desiredPoint = [ 2.5724 , 3.765 ]*freqScale;
% desiredPoint = [ 3.1724 , 2.965 ]*freqScale;
hold on
h = plot(desiredPoint(1),desiredPoint(2),'bd');
set(h,'MarkerFaceColor','r');
set(h,'MarkerSize',highlightMarkerSize)
hold off

%% Computing x-axis scaling in Frequency domain
nearestGridPoint = [round(desiredPoint(1))  round(desiredPoint(2))];
hold on
h = plot(nearestGridPoint(1),nearestGridPoint(2),'bo');
set(h,'MarkerFaceColor','b');
set(h,'MarkerSize',highlightMarkerSize)
hold off

if (nearestGridPoint(1) ~= 0)
    alpha1  =  desiredPoint(1) / nearestGridPoint(1);  % x axis
else
    alpha1 =  desiredPoint(1);
    nearestGridPoint(1) = 1 * sign(desiredPoint(1));
end


figure,
hold on
axis equal
if (alpha1 ~= 0)
     axis ([-N/2*alpha1 N/2*alpha1 -N/2 N/2])
end
for j  = 1: N+1                             % Marking X-Y grid lines
    h =plot(gridSpacing(:)*alpha1,gridSpacing(j)*ones(1,N+1),'b');
    h =plot(gridSpacing(j)*alpha1*ones(1,N+1),gridSpacing(:),'b');
end
for j  = 1: N+1                             % Plotting grid points
    pause(.03);
    h =plot(gridSpacing(j)*alpha1,gridSpacing(:),'bo');
    set(h,'MarkerFaceColor','g');
    set(h,'MarkerSize',gridMarkerSize)
end
h = plot(0,0,'bo');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize)

str = strcat ('X-axis 1D FrFT with scale, \alpha_1=', num2str(alpha1));
title(str)
xlabel ('u \rightarrow')
ylabel ('y \rightarrow')
hold off
 
hold on
h = plot(desiredPoint(1),desiredPoint(2),'bd');
set(h,'MarkerFaceColor','r');
set(h,'MarkerSize',highlightMarkerSize)
hold off

hold on
h = plot(nearestGridPoint(1)*alpha1,nearestGridPoint(2),'bo');
set(h,'MarkerFaceColor','b');
set(h,'MarkerSize',highlightMarkerSize)
hold off


%% Computing y-axis scaling in Frequency domain
if (nearestGridPoint(2) ~= 0)
    alpha2  =  desiredPoint(2) / nearestGridPoint(2);  % y axis
else
    alpha2 =  desiredPoint(2);
    nearestGridPoint(2) = 1 * sign(desiredPoint(2));
end

figure,
hold on
axis equal
if (alpha1 ~= 0 && alpha2 ~= 0)
     axis ([-N/2*alpha1 N/2*alpha1 -N/2*alpha2 N/2*alpha2])
end
for j  = 1: N+1                             % Marking X-Y grid lines
    h =plot(gridSpacing(:)*alpha1,gridSpacing(j)*ones(1,N+1)*alpha2,'b');
    h =plot(gridSpacing(j)*alpha1*ones(1,N+1),gridSpacing(:)*alpha2,'b');
end
for j  = 1: N+1                             % Plotting grid points
    pause(.03);
    h =plot(gridSpacing(j)*alpha1,gridSpacing(:)*alpha2,'bo');
    set(h,'MarkerFaceColor','g');
    set(h,'MarkerSize',gridMarkerSize)
end
h = plot(0,0,'bo');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize)

str = strcat ('Y-axis 1D FrFT with scale, \alpha_2=', num2str(alpha2));
title(str)
xlabel ('u \rightarrow')
ylabel ('v \rightarrow')
hold off

hold on
h = plot(desiredPoint(1),desiredPoint(2),'bd');
set(h,'MarkerFaceColor','r');
set(h,'MarkerSize',highlightMarkerSize)
hold off