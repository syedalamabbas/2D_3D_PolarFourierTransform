function [ DFT_Point ] = SinglePrecDirectSinglePoint_3DDFT( inputVolume, desiredPoint )
%COMPUTEFORSINGLEPOINT Summary of this function goes here
%   Computing the solution for a single point using row and column
%   differently scaled 1D FrFT
I = inputVolume;
[sizeX, ~, ~] =  size(I);             % Assume its a cube

N = sizeX -1;      % N is even

depthData = zeros(1, N+1,'single');

u = desiredPoint(1);
v = desiredPoint(2);
w = desiredPoint(3);

lineDataInner = zeros(1, N+1);
lineDataOuter = zeros(1, N+1);

n = -N/2: N/2;
Map_x = single(exp(-1i*2*pi*u*n/(N+1)));
Map_y = single(exp(-1i*2*pi*v*n/(N+1)));
Map_z = single(exp(-1i*2*pi*w*n/(N+1)));

% %% Fully direct computations - I
% for l = 1:N+1            % x
%     for m = 1: N+1       % y
%         depthData(:) = I (l,m,:);
%         lineDataInner(m) = sum(depthData .* Map_z);
%     end
%     lineDataOuter(l) =  sum(lineDataInner .* Map_x);
% end
% DFT_Point  = sum( lineDataOuter .* Map_y);


%% Fully direct computations - II
DFT_Point = single(0);
for l = 1:N+1            % x index
    for m = 1: N+1       % y index
        for k =1:N+1     % z index
            DFT_Point = DFT_Point + I(l,m,k)*single(exp(-1i*2*pi*(u*n(l)+v*n(m)+w*n(k))/(N+1)));   % Exactly according to the formula N^3 Complexity
        end
    end
end
end
