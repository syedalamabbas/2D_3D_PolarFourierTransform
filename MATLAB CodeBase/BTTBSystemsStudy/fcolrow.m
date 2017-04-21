%% A.9
function [ t ] = fcolrow( m,n, fchoice )
%FCOLROW generates the first columns and first rows of each block of the
%BTTB matrices T_{mn}. This program can also be used for nonsymmetric BTTB
% matrices
% m : the number of blocks;
% n : the size of each block
% fchoice : the choice of generating function;
% t : (2n)-by-(2m) matrix consists of first column and first rows of the
% blocks of the BTTB matric. Each column of t contains the first column and
% first row of a Toeplitz block. 

m1 = 2*m; n1 = 2*n; 
t = zeros (m1, n1);

for i = 0 :m-1
    for j= 0:n-1
        t(i+1,j+1)= kern(-i,-j, fchoice);
    end
end

for i = m+1:m1-1
    for j = 0:n-1
        t(i+1,j+1)= kern(m1-i, -j, fchoice);
    end
end

for i =0:m-1
    for j = n+1:n1-1
        t(i+1,j+1)= kern(-i, n1-j, fchoice);
    end
end

for i = m+1:m1-1
    for j= n+1:n1-1
        t(i+1,j+1) = kern(m1-i, n1 -j, fchoice);
    end
end

t = t.';

column1 = t(:,1)
figure, imagesc(toeplitz(column1(1:n1/2),column1(n1/2+1:n1)))
title('Toeplitz First column and first row in 1st column')
end

%% A.10

function y = kern(k,j,fchoice)
% kern computes the (j,1) entry of the (k,1) block of the BTTB matrix T_{mn}
% fchoice : the choice of the generating function: y : the output entry
j = abs(j); k = abs(k);  % for double symmetric BTTB matrix

if fchoice == 1         % Sequence (i)
    y = 1 ./ ( (j+1)*(k+1)^(1.1+0.1*j));
elseif fchoice == 2,    % Sequence (ii)
    y = 1 ./ ((j+1)^1.1* (k+1)^ (1.1+1.0*j));
elseif fchoice == 3,    % Sequence (iii)
    y = 1 ./ ((j+1)^1.1 + (k+1)^1.1);
elseif fchoice == 4,    % Sequence (iv)
    y = 1./ ((j+1)^2.1 + (k+1)^2.1);
end
end
