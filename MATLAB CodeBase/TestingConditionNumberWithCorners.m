clear
clc
close all
CondWithoutCorners = [];
CondWithCorners = [];
CondExtraPolar = [];
Condw_CornersExtraPolar = [];
CondPolar_Cartesian = [];
CondExtraPolar_Cartesian = [];
CondCartesian = [];
show = 1; % Turn it on to see the grid points but fix the ImageSize

ImageSizeRange = 24; 4:4:24;
count = 0;
for N = ImageSizeRange  % For different N
    M = 2*(N+2);
    
    %% Pure Cartesian
    if(show)
        figure,
        hold on
    end
    CartesianTransformMatrix = CreateMatrixFromGrid( 'C', N ,M,show );
    CondCartesian = [CondCartesian; cond(CartesianTransformMatrix )];
    if(show)
        xlabel('x')
        ylabel('y')
        axis equal
        axis ([-N/2 N/2 -N/2 N/2])
        hold off
    end
    
    %% Pure Polar x1
    if(show)
        figure,
        hold on
    end
    PurePolarTransFormMatrix =  CreateMatrixFromGrid( 'P', N ,M,show );
    CondWithoutCorners = [ CondWithoutCorners ; cond(PurePolarTransFormMatrix )];
    if(show)
        xlabel('x')
        ylabel('y')
        axis equal
        axis ([-N/2 N/2 -N/2 N/2])
        hold off
    end
    
    %% Polar + Corner x1
    if(show)
        figure,
        hold on
    end
    Polar_w_CornersTransFormMatrix =  CreateMatrixFromGrid( 'P_C', N ,M ,show);
    CondWithCorners = [ CondWithCorners; cond(Polar_w_CornersTransFormMatrix )];
    if(show)
        xlabel('x')
        ylabel('y')
        axis equal
        axis ([-N/2 N/2 -N/2 N/2])
        hold off
    end
    
    %% Pure Polar  x2
    if(show)
        figure,
        hold on
    end
    Polar_ExtraTransFormMatrix =  CreateMatrixFromGrid( 'E_P', N ,M,show );
    CondExtraPolar = [ CondExtraPolar; cond(Polar_ExtraTransFormMatrix )];
    if(show)
        xlabel('x')
        ylabel('y')
        axis equal
        axis ([-N/2 N/2 -N/2 N/2])
        hold off
    end
    
    %%  Polar + Corner  x2
    if(show)
        figure,
        hold on
    end
    Polar_w_CornersExtraTransFormMatrix =  CreateMatrixFromGrid( 'E_PC', N ,M, show);
    Condw_CornersExtraPolar = [ Condw_CornersExtraPolar; cond(Polar_w_CornersExtraTransFormMatrix )];
    if(show)
        xlabel('x')
        ylabel('y')
        axis equal
        axis ([-N/2 N/2 -N/2 N/2])
        hold off
    end
    
    %%  Polar + hybrid Corner  x1
    if(show)
        figure,
        hold on
    end
    Polar_CartesianTransFormMatrix =  CreateMatrixFromGrid( 'P_CartC', N ,M,show );
    CondPolar_Cartesian = [ CondPolar_Cartesian; cond(Polar_CartesianTransFormMatrix )];
    if(show)
        xlabel('x')
        ylabel('y')
        axis equal
        axis ([-N/2 N/2 -N/2 N/2])
        hold off
    end
    
    %%  Polar + hybrid Corner  x2
    if(show)
        figure,
        hold on
    end
    Polar_CartesianExtraTransFormMatrix =  CreateMatrixFromGrid( 'EP_CartC', N ,M,show);
    CondExtraPolar_Cartesian = [ CondExtraPolar_Cartesian; cond(Polar_CartesianExtraTransFormMatrix )];
    if(show)
        xlabel('x')
        ylabel('y')
        axis equal
        axis ([-N/2 N/2 -N/2 N/2])
        hold off
    end
    count = count +1 ;
    fprintf('\n Just finished iteration #%d',count);
end

figure,
semilogy (ImageSizeRange, CondCartesian, ImageSizeRange, CondWithoutCorners, ImageSizeRange, CondWithCorners, ImageSizeRange, CondExtraPolar,ImageSizeRange, Condw_CornersExtraPolar,ImageSizeRange, CondPolar_Cartesian,ImageSizeRange, CondExtraPolar_Cartesian, 'LineWidth', 1.7)
xlabel('N [Input array size is (N+1) \times (N+1)]');
ylabel('Condition-Number');
legend('Cartesian x1','Without Corners x1 Radially','With Corners x1 Radially', 'Without Corners x2 Radially','With Corners x2 Radially','With Corners x1 Radially Hybrid' , 'With Corners x2 Radially Hybrid' )
grid on
