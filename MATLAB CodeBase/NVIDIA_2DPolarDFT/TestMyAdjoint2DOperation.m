clear ;
clc;
close all

% %% PP  grid adjoint operation
% % <y,Ax> = <x,adj(A)y>'
% N=32; 
% X=randn(N)+sqrt(-1)*randn(N); Y=randn(2*N)+sqrt(-1)*randn(2*N);
% AX=PPFFT(X);
% AtY=APPFFT(Y);
% disp(abs( sum(sum(Y'.*AX)) -conj(sum(sum(X'.*AtY)))));
% 


%% MY Butterfly Grid Adjoint Operation
N = 88;   % always even
M = N+2;  % Always even again
X =  randn (N+1 );                            % Arbitrary Image Cartesian grid, this is a real image 
Y =  randn(M, N+1) + 1i* randn(M,  N+1);       % arbitrary Solution Polar grid , this must be complex values 
Ax = Compute2DPolarDFT (X, M);
AHy = Adjoint2DPolarDFT (Y);

prod1 = Y'.*Ax.';
prod2 = X'.*AHy;
disp(abs( sum(sum(prod1)) -conj(sum(sum(prod2)))));


%% This is for corners only
% 
% % checking using symbolic computations
 N = 66
 M = N+2 
% P = 1;       % Number of corner points derived from 1
% X = sym('X', [N+1,N+1])
% assume (X, 'real')
% % syms X1_1 X1_2 X1_3 X1_4 X1_5 X1_6 X1_7 X2_1 X2_2 X2_3 X2_4 X2_5 X2_6 X2_7 X3_1 X3_2 X3_3 X3_4 X3_5 X3_6 X3_7 X4_1 X4_2 X4_3 X4_4 X4_5 X4_6 X4_7 X5_1 X5_2 X5_3 X5_4 X5_5 X5_6 X5_7 X6_1 X6_2 X6_3 X6_4 X6_5 X6_6 X6_7 X7_1 X7_2 X7_3 X7_4 X7_5 X7_6 X7_7
% [Ax,~,C] = Compute2DPolarCornersDFT (X, M);
% % Ax = eval(Ax);
% 
% 
% % Forging matrix Ax = y where y is the stacked all polar corner points 
% columnStackedImage = reshape(X,[1,(N+1)^2]);
% zerosCell = cell(1, (N+1)^2-1 );
% for h = 1: (N+1)^2-1
%     zerosCell{h} = 0;
% end
% onesCell = cell(1,(N+1)^2 );
% for h = 1: (N+1)^2
%     onesCell{h} = 1;
% end
% indexes = 1:(N+1)^2;
% 
% syms MatA Temp  % the final A how it should look like
% counter = 1;
% for k = 1: P 
%     for c = 1:C
%         for z = 1: (N+1)^2                                 
%             returnVal = subs(Ax(k,c), columnStackedImage (indexes (indexes ~= z)), zerosCell ) ;
%             if(~isempty(returnVal))
%                  Temp(z) = returnVal;
%             else
%                 Temp(z) = 0;
%             end
%         end
%         InnerRow1 = (subs(Temp,columnStackedImage, onesCell));
%         InnerRow1  = eval(InnerRow1 );
%         fprintf('Just finished iteration #%d\n', counter);
%         MatA (counter,1: (N+1)^2 ) =InnerRow1; 
%         counter = counter+1;
%     end
% end
% 
% 
% Y1 = sym('Y', [P ,C])
% 
% [ AHy ] = Adjoint2DPolarCornersDFT( Y1, N,M );
% 
% 
% % Forging matrix AHy = x where y is the stacked all polar corner points and
% % X is the resultant Adjoint image
% columnStackedPolarCorners = reshape(Y1,[1,P *C]);
% zerosCell = cell(1, P *C-1 );
% for h = 1: P *C-1
%     zerosCell{h} = 0;
% end
% onesCell = cell(1,P *C );
% for h = 1: P *C
%     onesCell{h} = 1;
% end
% indexes = 1:P *C;
% 
% syms MatAH Temp  % the final A how it should look like
% counter = 1;
% for k = 1: N+1
%     for c = 1:N+1
%         for z = 1: P *C                                 
%             returnVal = subs(AHy(c,k), columnStackedPolarCorners (indexes (indexes ~= z)), zerosCell ) ; % Getting the image as columnwise stack
%             if(~isempty(returnVal))
%                  Temp(z) = returnVal;
%             else
%                 Temp(z) = 0;
%             end
%         end
%         InnerRow1 = (subs(Temp,columnStackedPolarCorners, onesCell));
%         InnerRow1  = eval(InnerRow1 );
%         fprintf('Just finished iteration #%d\n', counter);
%         MatAH (counter,1: P*C  ) =InnerRow1; 
%         counter = counter+1;
%     end
% end  
% 
% MatA = eval(MatA); 
% MatAH = eval(MatAH);
% error = MatA - MatAH'                        % This should be zero for the adjoint to work
% MatA(1,:).' - conj(MatAH(:,1) )
% MatA(2,:).' - conj(MatAH(:,2) )

X =  randn (N+1 );                             % Arbitrary Image Cartesian grid, this is a real image
% X =  zeros (N+1 );   
% X(1) = 1;         % Checking with impulse
[Ax,~,C] = Compute2DPolarCornersDFT2(X, M);
figure, imagesc(X)

% returnImage = Adjoint2DPolarCornersDFT (Ax, N,M);
% figure, imagesc(real(returnImage),[0 255])
%  
Y =  randn(4, C) + 1i* randn(4, C);       % arbitrary values Polar grid corners, this must be complex values 
AHy = Adjoint2DPolarCornersDFT2(Y, N,M);

prod1 = Y'.*Ax.';  
prod2 = X'.*AHy;
% disp(eval(abs( sum(sum(prod1)) -conj(sum(sum(prod2)))))); % Use this for symbolic
disp((abs( sum(sum(prod1)) -conj(sum(sum(prod2)))))); 

 

%% What is the system of corners

X =  zeros (N+1 );   
X(1) = 1; 
VectorFirstCol2D = Adjoint2DPolarCornersDFT2(Compute2DPolarCornersDFT2(X, M), N,M);

figure, imagesc(real(VectorFirstCol2D))





