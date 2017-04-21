
clear
clc
close all

%% Preparing the image
Size = 251;                             % Odd
% I = double(imresize(imread('lena.tif'), [Size,Size]));
I = double(imresize(imread('cameraman.tif'), [Size,Size]));

figure , imshow(I, [0 255]);
title ('Original Image')
N = Size -1;                       % N is even


%% Padding the image
NewN_2 = floor(N/2/(cosd(45)));
NewSize = 2*NewN_2+1;

I_New =  padarray (I, [NewN_2-N/2, NewN_2-N/2], 'both' );
figure , imshow(I_New, [0 255]);
title ('Padded Image')

NewSize = Size;
I_New = I;

%% Taking the Polar DFT
noOfAngles = (NewSize +1);
M = noOfAngles;
Polar_Grid = VectorizedCompute2DPolarDFT( I_New,  M ) ;

%% Plotting it
deltaTheta = 180/M;     % Angular sampling rate
angles = [0:deltaTheta:180-deltaTheta];  
UniformSpacing1 = [-N/2 N/2];
UniformSpacing2 = [-NewN_2 NewN_2];
figure, imagesc(angles,UniformSpacing2,log(abs(Polar_Grid')))
title('F_{\theta} (X)');
xlabel('\theta (degrees)');
ylabel('X'); 
colormap(hot);
colorbar

% %% Nulling the unwanted frequencies 
% for k = 1:length(angles)
%     angle = angles(k);
%     if (  angle <= 45 || angle > 135 )
%         newGridSpacing = UniformSpacing1 ./ cosd(angle);       % BH
%         if(angle > 135)
%             halfSpacing = 1:1:newGridSpacing(1);
%         else
%             halfSpacing = 1:1:newGridSpacing(end);
%         end
%     else if ( angle > 45 || angle <= 135 )
%             newGridSpacing = UniformSpacing1 ./ sind(angle);    %BV
%             halfSpacing = 1:1:newGridSpacing(end);
%         end
%     end
%     
%     %  newGridSpacing = [-fliplr(halfSpacing) 0 halfSpacing]; <=   UniformSpacing2 
%     
%     CurrentRadialLine = Polar_Grid(k,:);                                                  % Retrieving 
%     indexesToNull = [ (-[halfSpacing(end)+1:NewN_2]+ NewN_2+1)  ([halfSpacing(end)+1:NewN_2]+ NewN_2+1)  ] ;       % Processing
%     CurrentRadialLine(indexesToNull) = zeros(1,length(indexesToNull));                    % Updating
%     Polar_Grid(k,:) =  CurrentRadialLine ;                                                 % Replacing           
% end 
% 
% figure, imagesc(angles,UniformSpacing2,log(abs(Polar_Grid')))
% title('Nulled Frequencies outside the square Image, F_{\theta} (X)');
% xlabel('\theta (degrees)');
% ylabel('X'); 
% colormap(hot);
% colorbar


%% Computing the inverse with those nulled frequencies
[ ReconstImage ] = Inverse2DPolarDFT( Polar_Grid, N/2 );

ReconstImage = real(ReconstImage);
% R = [0 255];
% dR = diff( R );
% 
% ReconstImage =  ReconstImage - min( ReconstImage(:)); % set range of A between [0, inf)
% ReconstImage = ReconstImage  ./ max( ReconstImage(:)) ; % set range of A between [0, 1]
% ReconstImage=  ReconstImage .* dR ; % set range of A between [0, dRange]
% ReconstImage =  ReconstImage + R(1); % shift range of A to R

figure,     
 hold on
subplot(2,2,1) % first subplot    
axis ([-Size Size -Size Size]) 
imshow(I_New, [0 255]);
% imshow(I)
xlabel('x')
ylabel('y');  
title('Original Image')  
subplot(2,2,2) % second subplot
axis equal
imshow(real(ReconstImage), []);  
xlabel('x')  
ylabel('y');
title('Reconstructed Image fast exact solution')
disp('The difference norm between the exact and reconstructed is ')
disp(norm(real(ReconstImage)-I_New)); 
  
figure, imagesc((ReconstImage-I_New)) ; colorbar
title('Reconstruction Error as fixed pattern noise')



