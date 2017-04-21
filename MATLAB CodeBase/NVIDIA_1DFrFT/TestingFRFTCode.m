clear all
clc;
close all
%% This script is written purely for debugging purposes 
% We TEST various functions through direct computations 
% and through fast implemented version to mostly check the accuracy, follow
% along

%% Declaring few variables
N = 1351;   % always odd 
x = double(rand(1, N)');
%  x = zeros(1,N)';     % This is the delta function to be used
%  x(3) = 1;
alpha = double(.9569);

 
i=sqrt(-1);                 % initializing incase it hold another value

%% Direct FRFT
y1=zeros(N,1);
for n=-N/2:1:N/2-1,
    for k=0:1:N-1,
        y1(n+1+N/2)=y1(n+1+N/2)+x(k+1)*exp(-i*2*pi*k*n*alpha/N);
    end;
end;
 
  
%% The convolution via FFT that is centered only in Frequency Domain
n=[0:1:N-1, -N:1:-1]';
Factor1=exp(-i*pi*alpha*n.^2/N);
Factor2=exp(i*pi*(0:1:N-1)'*N*alpha/N);

x_tilde=x.*Factor2;
x_tilde=[x_tilde; zeros(N,1)];
x_tilde=x_tilde.*Factor1;    

XX=fft(x_tilde);
YY=fft(conj(Factor1));
y4=ifft(XX.*YY);
y4=y4.*Factor1;
y4=y4(1:N);

disp ('Is it working General Centered ???')
disp(max(abs(y1-y4))); 


%% Direct FRFT fully Centered in both spatial and frequency domains
m =1;
Centeredy1=zeros(N,1);
for n=-(N-1)/2:(N-1)/2   
    z = 1;
    for k=-(N-1)/2:1:(N-1)/2
        Centeredy1(m)=Centeredy1(m)+x(z)*exp(-i*2*pi*k*n*alpha/N);
        z = z+1;
    end;
    m = m+1;
end; 

% Computing inverse
m =1;
ICenteredy1=zeros(N,1);
for n=-(N-1)/2:(N-1)/2   
    z = 1;
    for k=-(N-1)/2:1:(N-1)/2
        ICenteredy1(m)=ICenteredy1(m)+Centeredy1(z)*exp(i*2*pi*k*n*alpha/N);
        z = z+1;
    end;
    m = m+1;
end;

% Convolution via FFT Centered Function
n=[0:1:N-1, -N:1:-1]';

Factor1=exp(-i*pi*alpha*n.^2/N) ;

M = (0:1:N-1);
M = M - (N-1)/2;
Factor2=exp(i*pi*M'*(N-1)*alpha/N);

k = (0:1:N-1)';
Factor3 = exp(i*pi*alpha*(N-1)*k/N);

x_tilde=x.*Factor2;
x_tilde=[x_tilde; zeros(N,1)];
x_tilde=x_tilde.*Factor1;    

XX=fft(x_tilde);
YY=fft(conj(Factor1));
Newy4=ifft(XX.*YY);
Newy4=Newy4.*Factor1;
Newy4=(Factor3.*Newy4(1:N));   

disp ('Is it working Centered explicit computations???')
disp(max(abs(Centeredy1-Newy4)));


%--------------------------- My Functions test---------------------------%

%% Check the direct implementation 
Centeredy11 = Direct1DFrFT(x' ,alpha );                  % Direct Check
disp ('Is it working Centered my function Direct1DFrFt  ???')
disp(max(abs(Centeredy1- reshape(Centeredy11, [N,1]))));




%% Check my function
NewFFTCentered = FrFT_Centered(x,alpha);      % 9 x 1
disp ('Come on, is it finally working Centered  my function FrFT_Centered  ???')
disp(max(abs(Centeredy1- NewFFTCentered)));


NewFFTCentered = FrFT_Centered(x',alpha);      % 1 x 9
disp ('Come on, is it finally working Centered  my function FrFT_Centered  ???')
disp(max(abs(Centeredy11- NewFFTCentered)));


% Computed from single points using FrFTSingleCentered code
desiredIndexes = [-(N-1)/2 : (N-1)/2];
lineData = zeros(N,1);
for i = 1:N
    lineData(i) = FrFTCenteredSingle( x, alpha, desiredIndexes(i) );
end

disp ('Is single point  my function FrFTCenteredSingle working  ???')
disp(max(abs(Centeredy1- lineData)));

%% Checking Vectorized
vectSol  = VectorizedFrFT_Centered(x,alpha);
ctime = timeit(@()VectorizedFrFT_Centered(x,alpha));     % Measure CPU time
disp(['Execution time on CPU (Vectorized) (in seconds)= ',num2str(ctime)]);
disp ('Maximum absolute error forward FrFT:')
disp(max(abs(Centeredy1- vectSol)));  
  
%  %% Checking the accuracy of the inverse computations
% InverseReconstruction = VectorizedFrFT_Centered(vectSol,-alpha);
% 
%  
% figure, plot ( 1:N, x, 1:N, real(ICenteredy1)/ N , '+g', 1:N, real(InverseReconstruction)/N, 'LineWidth', 2.1)
% legend('original', 'reconstruction direct',  'reconstruction myFrFT')
% grid on
% 
% figure, plot ( 1:N,abs( x - real(ICenteredy1)/ N), '+g',  1:N,abs( x - ifft(fft(x))) , 1:N,abs( x - real(InverseReconstruction)/N), 'LineWidth', 2.1)
% legend('Error using direct computation FrFT' , 'error using FFT', 'error using my code FrFT')
% grid on
% disp ('Maximum absolute error Inverse FrFT same alpha and vector')
% disp(max(abs(ICenteredy1- InverseReconstruction)));
%    
%  
% %% Checking GPU ! Precomputations again on the GPU
% gx = gpuArray(double(x));             % Move to GPU
% [sizeX,sizeY] =  size(x);
% N = sizeX -1;      % N is even
% M = N+1;                                     % Number of angles
% L = 1;
% 
% n=[0:1:N, -N-1:1:-1]';
% J = (0:1:N)';
% K = (0:1:N)';
% Z = zeros(N+1,1,'like',gx);
% % Z = zeros(N+1,N+1, 'like',I);
% J = J - N/2;
% 
% E_n =  exp(-1i*pi*alpha*n.^2/(N+1));                  % Sequence as defined in the paper
% PremultiplicationFactor = exp(1i*pi*J *N*alpha/(N+1));
% PostmultiplicationFactor = exp(1i*pi*alpha*N*K/(N+1));
% 
% gpuSol = GPUVectorizedFrFT_Centered(gx, Z, E_n, PremultiplicationFactor, PostmultiplicationFactor);
% gtime = gputimeit(@()GPUVectorizedFrFT_Centered(gx, Z, E_n, PremultiplicationFactor, PostmultiplicationFactor));% Measure GPU time
% disp(['Execution time on GPU (Vectorized) (in seconds) = ',num2str(gtime)]);
% disp ('Maximum absolute error')
% disp(max(abs(Centeredy1 - gather(gpuSol))));

 






