function [ PolarGrid ] = ComputeSymbolic2DPolarDFT( Image,  M , alphas, betas, one_alpha)
%COMPUTESYMBOLIC2DPOLARDFT Summary of this function goes here
%   Detailed explanation goes here
[Size,~] = size(Image);
spacing = -(Size-1)/2 : (Size-1)/2;
IndexBetas = spacing' * spacing ;

N = Size -1;
L = (M-2)/4;                                   % Number of levels
if(rem(M-2,4) ~= 0)
    L = ceil (L);                % use + 1 to compute for 45
end

for l=1:1:L  % For each level
    angle = l*180/M;
    alpha_l =  alphas(l);
    beta_l  =  betas(l);
    
    % X axis Scaling
    for p = 1: Size      % For row
        for q = 1: Size   % For columns
            FrFT_XAxis_Uniform(p,q) = SymbolicFrFT_Uniform( Image(p,:), alpha_l, spacing(q), Size );
        end
    end
    Line_1 = sum ( FrFT_XAxis_Uniform .* exp(-2*1i*pi*IndexBetas* beta_l/ Size)  );
    Line_1Conj = sum ( FrFT_XAxis_Uniform .* exp(2*1i*pi*IndexBetas* beta_l/ Size)  );
    
    PolarGrid ( 1+l, 1:Size)= Line_1;
    PolarGrid (M+1-l, 1:Size)= fliplr( Line_1Conj);
    
    if (angle == 45)   % If angle is 45 we have already computed it in X butterfly grid
        continue;
    end
    
    % Y axis Scaling
    for p = 1: Size      % For row
        for q = 1: Size   % For columns
            FrFT_YAxis_Uniform(p,q) = SymbolicFrFT_Uniform( Image(:,p).', alpha_l, spacing(q), Size );
        end
    end
    Line_2 = sum ( FrFT_YAxis_Uniform .* exp(-2*1i*pi*IndexBetas* beta_l/ Size)  );
    Line_2Conj = sum ( FrFT_YAxis_Uniform .* exp(2*1i*pi*IndexBetas* beta_l/ Size)  );
    
    PolarGrid (M/2+1+l, 1:Size)=   Line_2Conj ;
    PolarGrid (M/2+1-l, 1:Size)=   Line_2;
                                                          
    if (l == 1)
%         Temp = .';
        NintyLine = FrFT_XAxis_Uniform( 1:Size,N/2+1).'; %Temp(N/2+1,1:Size);                 % Getting the column of 90 degrees angle
        ZeroLine = FrFT_YAxis_Uniform( 1:Size,N/2+1).';   % Getting the row of 0 degrees angle, it is rotated since it has the values that are scaled in one axis
        
        for q = 1: Size   % For columns
            DFT2_ZeroLine (q) =  SymbolicFrFT_Uniform( ZeroLine, one_alpha, spacing(q), Size );
            DFT2_NinetyLine(q) = SymbolicFrFT_Uniform( NintyLine, one_alpha, spacing(q), Size );
        end
        PolarGrid (M/2+1, 1:Size) =   DFT2_NinetyLine;
        PolarGrid (1, 1:Size) = DFT2_ZeroLine;
    end
end
end