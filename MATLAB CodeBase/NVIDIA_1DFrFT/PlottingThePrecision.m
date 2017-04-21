close all
clear 
clc


%% Some initializations  Declaring few variables
i=sqrt(-1);                 % initializing incase it holds another value

%% Iterating through the loop
H = 300;
factor = 5;
N_0 = 10;

Max_Abs_Error = zeros(H,1);
Max_Abs_Error_Split = zeros(H,1);
for h = 1: H
    
    %% Initializations
    N = N_0+h*factor;          %% arithmetic progression
    x = double(10*rand(1, N)');
    alpha = double(.9569);
%     alpha = double(rand(1,1));
    
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
    end
    
    %% Checking Vectorized
    vectSol  = VectorizedFrFT_Centered(x,alpha);
    disp(['Iteration:',num2str(h)]);
    disp ('Maximum absolute error forward FrFT:');
    Max_Abs_Error(h) =  max(abs(Centeredy1- vectSol));
    disp(Max_Abs_Error(h));
    
    
    %% Checking Vectorized with split
    vectSol2  = SplitVectorizedFrFT_Centered(x,alpha);
    disp(['Iteration:',num2str(h)]);
    disp ('Maximum absolute error forward FrFT:');
    Max_Abs_Error_Split(h) =  max(abs(Centeredy1- vectSol2));
    disp(Max_Abs_Error_Split(h));
end
 

%% Plotting final results
figure, 
semilogy(N_0+(1:H)*factor,  Max_Abs_Error,'LineWidth', 2.3)
% semilogy(N_0+(1:H)*factor,  Max_Abs_Error,'r', N_0+(1:H)*factor, Max_Abs_Error_Split,'LineWidth', 2.3)
xlabel('Vector length, N\rightarrow')
ylabel('Max Abs Error')
legend('1D FrFT using fast convolution')
% legend('Vectorized full length', 'Vectorized with high padding')
grid on