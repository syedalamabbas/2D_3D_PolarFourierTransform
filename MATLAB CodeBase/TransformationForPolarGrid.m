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
 title('\bf{ Polar grid}') 
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

%% Is there a reverse solution possible ???

%% Invertible Polar Grid
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
    set(h,'MarkerSize',gridMarkerSize-2)
end
h = plot(0,0,'bs');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize)

title('\bf{Stably Invertible Polar grid}') 
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
hold off

M = N + 2;                 % Total number of polar lines 
deltaTheta = 180/M;     % Angular sampling rate
angles = [0:deltaTheta:180];  

newGridSpacing = 0;                    % This is a change
hold on  
for angle = angles
    pause(delay)
    
    if (  angle <= 45 || angle > 135 )
        newGridSpacing = gridSpacing ./ cos(angle*pi/180);       % BH
        if(angle > 135)
            halfSpacing = 1:1:newGridSpacing(1);
        else
            halfSpacing = 1:1:newGridSpacing(end);
        end
    else if ( angle > 45 || angle <= 135 )
            newGridSpacing = gridSpacing ./ sin(angle*pi/180);    %BV
            halfSpacing = 1:1:newGridSpacing(end);
        end
    end
    
    newGridSpacing = [-fliplr(halfSpacing) 0 halfSpacing];
    
    h1= plot(newGridSpacing*cosd(angle),newGridSpacing*sind(angle),'ko');
    set(h1,'MarkerFaceColor','r');
    set(h1,'MarkerSize',gridMarkerSize)
    
    halfPoints = halfSpacing(N/2+1:end);    % These are the corner points
    if(~isempty(halfPoints))
        halfPoints = [-halfPoints ; halfPoints];
        for point = halfPoints
        h1= plot(point *cosd(angle),point *sind(angle),'ko');
        set(h1,'MarkerFaceColor','g');
        set(h1,'MarkerSize',gridMarkerSize+1)
        end
    end
    
end 
hold off    
   

%% Second Invertible Polar Grid -My type
figure,axis equal
hold on
axis equal
for j  = 1: N+1                             % Marking X-Y grid lines
    h =plot(gridSpacing(:),gridSpacing(j)*ones(1,N+1),'b');
    h =plot(gridSpacing(j)*ones(1,N+1),gridSpacing(:),'b');
end
for j  = 1: N+1                             % Plotting grid points
    pause(delay);
    h =plot(gridSpacing(j),gridSpacing(:),'bs');
    set(h,'MarkerFaceColor','g');
    set(h,'MarkerSize',gridMarkerSize-2)
end
h = plot(0,0,'bs');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize)

title('\bf{Stably Invertible Polar grid II}') 
xlabel ('\bf{x \rightarrow}')
ylabel ('\bf{y \rightarrow}')
hold off

newGridSpacing = gridSpacing ./cosd(45);
 halfSpacing = 1:1:newGridSpacing(end);
 newGridSpacing = [-fliplr(halfSpacing) 0 halfSpacing];

 new__N = halfSpacing(end);
 axis ([-new__N new__N -new__N new__N])
hold on 
for angle = angles 
    pause(delay) 
    h1= plot(newGridSpacing*cosd(angle),newGridSpacing*sind(angle),'ko');
    set(h1,'MarkerFaceColor','r');
    set(h1,'MarkerSize',gridMarkerSize) 
end
hold off   

%% First Scale sind(18)
alpha11 =  cosd(deltaTheta); 
% alpha11 = sind(deltaTheta);           % Switching the scales to see geometry
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
% beta11 = sind(3*deltaTheta);                            % Switching the scales to see geometry
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
    set(h,'MarkerSize',gridMarkerSize+4)
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
 
%% Third Squeez sind(3*18) 
alpha13 = alpha11; 
figure,axis equal  
hold on
axis equal
axis ([-N/2 N/2 -N/2 N/2])
for j  = 1: N+1                             % Marking X-Y grid lines
    h =plot(gridSpacing(:),gridSpacing(j)*ones(1,N+1)*alpha13,'b');
    h =plot(gridSpacing(j)*ones(1,N+1),gridSpacing(:)*alpha13,'b');
end
for j  = 1: N+1                             % Plotting grid points
    pause(delay);
    h =plot(gridSpacing(j),gridSpacing(:)*alpha13,'bs');
    set(h,'MarkerFaceColor','g');
    set(h,'MarkerSize',gridMarkerSize)
end
h = plot(0,0,'bs');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize)

title('\bf{1^{st} level scaling, y-axis by cos(\Delta\theta)}')
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

%% Next scale in x direction
figure,axis equal
hold on
axis equal
% axis ([-N/2 N/2 -N/2 N/2])
for j  = 1: N+1                             % Marking X-Y grid lines
    h =plot(gridSpacing * newspacing (j) ,gridSpacing(j)*ones(1,N+1)*alpha13,'b');
    h =plot(gridSpacing(j)*ones(1,N+1)* newspacing(j),gridSpacing(:)*alpha13,'b');
end
for j  = 1: N+1                             % Plotting grid points
    pause(delay);
    h =plot(gridSpacing(:)* newspacing(j),gridSpacing(j)*alpha13,'bs');
    set(h,'MarkerFaceColor','g');
    set(h,'MarkerSize',gridMarkerSize)
end
h = plot(0,0,'bs');                 % DC value
set(h,'MarkerFaceColor','k'); 
set(h,'MarkerSize',gridMarkerSize)

title('\bf{1^{st} level scaling, y-axis by cos(\Delta\theta) then in x-axis by sin(\Delta\theta)}')
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
axis tight

 
%% Second Squeez sind(2*18)
alpha12 =  cosd(deltaTheta*2);
figure,axis equal
hold on
axis equal
axis ([-N/2 N/2 -N/2 N/2]) 
for j  = 1: N+1                             % Marking X-Y grid lines
    h =plot(gridSpacing(:)*alpha12,gridSpacing(j)*ones(1,N+1),'b');
    h =plot(gridSpacing(j)*ones(1,N+1)*alpha12,gridSpacing(:),'b');
end
for j  = 1: N+1                             % Plotting grid points
    pause(delay);
    h =plot(gridSpacing(j)*alpha12,gridSpacing(:),'bs');
    set(h,'MarkerFaceColor','g');
    set(h,'MarkerSize',gridMarkerSize)
end
h = plot(0,0,'bs');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize)

title('\bf{2^{nd} level scaling, x-axis by cos(2\Delta\theta)}')
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
beta12 = sind(2*deltaTheta);
newspacing = abs(gridSpacing) * beta12/2;

newspacing(N/2+1) = 1;          % Scaling reset at zero

figure,axis equal
hold on 
axis equal
% axis ([-N/2 N/2 -N/2 N/2])
for j  = 1: N+1                             % Marking X-Y grid lines
    h =plot(gridSpacing(:)*alpha12,gridSpacing(j)*ones(1,N+1)*newspacing(j),'b');
    h =plot(gridSpacing(j)*ones(1,N+1)*alpha12,gridSpacing(:)*newspacing(j),'b');
end
for j  = 1: N+1                             % Plotting grid points
    pause(delay+.02);
    h =plot(gridSpacing(j)*alpha12,gridSpacing(:)*newspacing(j),'bs');
    set(h,'MarkerFaceColor','g');
    set(h,'MarkerSize',gridMarkerSize)
end
h = plot(0,0,'bs');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize)

title('\bf{2^{nd} level scaling, x-axis by cos(2\Delta\theta) then y axis by sin(2\Delta\theta)}')
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
axis tight

%% Fourth Squeez sind(4*18) 
alpha14 =  alpha12; 
figure,axis equal
hold on
axis equal
axis ([-N/2 N/2 -N/2 N/2])
for j  = 1: N+1                             % Marking X-Y grid lines
    h =plot(gridSpacing(:),gridSpacing(j)*ones(1,N+1)*alpha14,'b');
    h =plot(gridSpacing(j)*ones(1,N+1),gridSpacing(:)*alpha14,'b');
end
for j  = 1: N+1                             % Plotting grid points
    pause(.03);
    h =plot(gridSpacing(j),gridSpacing(:)*alpha14,'bs');
    set(h,'MarkerFaceColor','g'); 
    set(h,'MarkerSize',gridMarkerSize)
end
h = plot(0,0,'bs');                 % DC value
set(h,'MarkerFaceColor','k');
set(h,'MarkerSize',gridMarkerSize)

title('\bf{2^{nd} level scaling, y-axis by cos(2\Delta\theta)}')
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

%% Next scale in x direction
figure,axis equal
hold on
axis equal
% axis ([-N/2 N/2 -N/2 N/2])
for j  = 1: N+1                             % Marking X-Y grid lines
    h =plot(gridSpacing * newspacing (j) ,gridSpacing(j)*ones(1,N+1)*alpha14,'b');
    h =plot(gridSpacing(j)*ones(1,N+1)* newspacing(j),gridSpacing(:)*alpha14,'b');
end
for j  = 1: N+1                             % Plotting grid points
    pause(delay);
    h =plot(gridSpacing(:)* newspacing(j),gridSpacing(j)*alpha14,'bs');
    set(h,'MarkerFaceColor','g');
    set(h,'MarkerSize',gridMarkerSize)
end
h = plot(0,0,'bs');                 % DC value
set(h,'MarkerFaceColor','k'); 
set(h,'MarkerSize',gridMarkerSize)

title('\bf{2^{nd} level scaling, y-axis by cos(2\Delta\theta) then in x-axis by sin(2\Delta\theta)}')
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
axis tight



% %% Plotting irrational numbers as period
% 
% z = [1:(N+1)+3];
% alpha =  cosd(18);
% figure, stem (z,sin(2*pi*z*alpha/(N+1)), 'Linewidth', 2.3)
% axis( [ 1 (N+1)+3  -1.5 1.5 ] )
% xlabel('\bf{ n \rightarrow}')
% ylabel(' \bf{  sin(2\pi n \alpha /(N+1)) \rightarrow }')
% grid on
% % legend('\bf{cos(\Delta \theta)}')
