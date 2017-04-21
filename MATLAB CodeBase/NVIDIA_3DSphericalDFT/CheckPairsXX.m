function IsPresent = CheckPairsXX(k,m, K,M )

%% Extra coding for XX ,YY , XZ and YZ blocks seperate conditional computations
% No of levels theta                           % These computed on fly
P = (K-2)/4;                                   % Number of levels
if(rem(K-2,4) ~= 0)
    P = ceil (P);                           % use + 1 to compute for 45 but skip for other dimension
end

% No of levels Phi
Q = (M-2)/4;                                   % Number of levels
if(rem(M-2,4) ~= 0)
    Q = ceil (Q);                           % use + 1 to compute for 45 but skip for other dimension
end

persistent XXpairs;
if(isempty(XXpairs))
    for q=1:Q  % For each level of Polar slices
        %% X oriented pair of Polar slices, (1+q, K+1-q)
        for p = 1:P
%             if(p == 1 && q == 1)              % First lines only
                XXpairs = [XXpairs;  p+1, q+1 ; ];
                XXpairs = [XXpairs;  p+1, M+1 - q;];
                XXpairs = [XXpairs;  K+1 - p, q+1;];
                XXpairs = [XXpairs;  K+1 - p, M+1 - q;];
%             end
        end
    end
end
[sizeX,~ ] = size(XXpairs);
IsPresent = 0;
for g=1:sizeX
    currentPairs = XXpairs(g,:);
    angleTheta_Index = currentPairs(1);
    anglePhi_Index = currentPairs(2);
    if(angleTheta_Index == k && anglePhi_Index == m)
        IsPresent = 1;
        break;
    end
end

end