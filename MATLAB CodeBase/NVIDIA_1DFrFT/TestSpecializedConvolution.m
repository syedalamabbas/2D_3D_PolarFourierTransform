
close all
clear
clc

H = 500;
factor = 5;
N_0 = 10;

Max_Abs_Error_Special = zeros(H,1);
Max_Abs_Error_fast = zeros(H,1);
for h = 1: H
     
    N = N_0+h*factor;
    f = complex(rand(N,1),rand(N,1));
    g = complex(rand(N,1),rand(N,1));
    
    directconv = conv(f,g);
    cconv = SpecialDealiasedConvolution( f, g );
    Max_Abs_Error_Special(h) = max(abs(directconv(1:N)-cconv));
    
    
    % These three operations constitute "fast convolution" or FFT based convolution
    f = [f; zeros(size(f))];
    g = [g; zeros(size(g))];
    FirstFFT=fft(f);
    SecondFFT=fft(g);
    fastconv=ifft( bsxfun(@times,FirstFFT,SecondFFT));             % Computing Convolution here
    Max_Abs_Error_fast(h) = max(abs(directconv(1:N)-fastconv(1:N)));
end

%% Simple Comparison plot
figure, plot(real(directconv), 'r')
hold on
plot(real(cconv), 'g')
plot(real(fastconv),'b')
hold off
xlabel('Data points')
ylabel('Amplitude')
legend('Direct','SpecialDealiased')
grid on

disp('Max Absolute error first half')
disp(max(abs(directconv(1:N)-cconv)))


%% Plotting final results
figure,
semilogy(N_0+(1:H)*factor,  Max_Abs_Error_Special,N_0+(1:H)*factor,  Max_Abs_Error_fast,'LineWidth', 2.3)
xlabel('Vector length, N\rightarrow')
ylabel('Max Abs Error')
legend('Specialized Direct', 'Simple Fast with padding')
grid on
