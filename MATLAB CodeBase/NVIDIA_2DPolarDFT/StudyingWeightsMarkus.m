clear 
clc


M = 128;   % Number of angles
N = 64;   % Number of points
TotalNoOfNodes = M*(N+1);

%% According to MArkus 2007 paper
W = pi*abs(-N/2:N/2)/(M*(N+1)^2);
W(N/2+1) = pi/(4*M*(N+1)^2);

Mat = ones(M,1)*W;
figure, imagesc(Mat')
colorbar
title(' MArkus 2007 paper')

%% Accordiung to code
Const_factor=M*(((N+1)/2.0)*((N+1)/2.0)+1.0/4.0);
W = abs(-N/2:N/2)/Const_factor;
W(N/2+1) = 1.0/4.0/Const_factor;

Mat = ones(M,1)*W;
figure, imagesc(Mat')
colorbar
title(' MArkus 2007 paper code ')

%% Amir
W =[];
W         = sqrt(abs(-N/2:N/2)/2)/(N+1);    % Fourier based preconditioner from Amir's 2001 paper Fast Slant Slack
W(N/2+1)  = sqrt(1/8)/(N+1);


Mat = ones(M,1)*W; 
figure, imagesc(Mat'.* Mat')
colorbar
title(' Amir 2007 paper code ')
