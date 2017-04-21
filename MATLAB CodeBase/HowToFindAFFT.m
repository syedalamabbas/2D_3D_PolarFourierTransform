 figure,
 axis([-5, 5, -5 , 5 ])
hold on
plot([-4:4],ones(1,9), '+r', 'LineWidth',1.7)  
plot([-3.5:1:3.5],2*ones(1,8), '+b', 'LineWidth',1.7)  
plot([-4:4],0*ones(1,9), '+g', 'LineWidth',1.7)    
plot([-4:4]*(3.5/4),-1*ones(1,9), 'ok', 'LineWidth',1.7)  
hold off
grid on
legend('Spatial Points','Desired FFT','Normal FFT', 'Fractional FFT, \alpha = 0.8750' )