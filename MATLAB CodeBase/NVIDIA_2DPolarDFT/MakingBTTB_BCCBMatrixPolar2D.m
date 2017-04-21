clear;
clc
close all


% %% checking that you can compute the column of TBBT without computing the product A^H A explicitly
% N = 4;  % always even
% M = N+2 ;  % always even
% 
% U = sym('u', [M*(N+1)^2  (N+1)^2]);
% V = sym('v', [M*(N+1)  M*(N+1)^2]);
% % x = sym('x', [(N+1)^2 (N+1)^2]);
% 
% 
% fullProd = U'*V'*V*U;
% 
% 
% x_valMat = zeros((N+1)^2,(N+1)^2);
% x_valMat(1,1) = 1; 
% 
% directComp = U'*V'*V*U*x_valMat;
% directComp(:,1) - fullProd(:,1)

N = 10;     % Always Even
M = 1*(N+2);   % Always Even


%% Creating the BTTB matrix From the first column only 
ImageFinal = zeros(N+1,N+1);
ImageFinal(1,1) = 1;            % Making it an impulse or a point spread function, 

W         = sqrt(abs(-N/2:N/2)/2)/(N+1);    % Fourier based preconditioner from Amir's paper
W(N/2+1)  =  sqrt(1/8)/(N+1);
W        = ones(M,1)*W ;

BTTB_ColumnsIn2DShape = Adjoint2DPolarDFT(W.*W.* Compute2DPolarDFT( ImageFinal,  M ));     % It is (N+1) x (N+1) size, it is a impluse response of the linear system 
%  BTTB_ColumnsIn2DShape = Adjoint2DPolarDFT(Compute2DPolarDFT( ImageFinal,  M ));     % Without preconditioner W

% % For symbolic  do this  ---->
% N = 4;
% mat_bttb = sym('a_', [N+1,N+1]);
% assume(mat_bttb, 'real')
% BTTB_ColumnsIn2DShape = mat_bttb;


BTTB_ColumnsIn2DShape = real(BTTB_ColumnsIn2DShape);
Indexes = 1:(N+1);
Toeplitz_Indexes = toeplitz(Indexes);   % Block Level

BTTB_Matrix = []; 
for l = 1:(N+1)
    ColumnIndexes = Toeplitz_Indexes(:,l);
    ToeplitzBlocks = [];
    for k = 1:(N+1)
        CurrentBlockIndex = ColumnIndexes(k);                      % This is the block Index
        ToeplitzBlock = toeplitz(BTTB_ColumnsIn2DShape(:,CurrentBlockIndex));
        ToeplitzBlocks = [ ToeplitzBlocks ; ToeplitzBlock];    % Append Column-wise T = [T1 ; T2; T3]
    end
    BTTB_Matrix  = [BTTB_Matrix, ToeplitzBlocks];              % Append Row-wise
end

% Regularized matrix
alpha = .5;
BTTB_Matrix = BTTB_Matrix + eye((N+1)^2)*alpha ;

BTTB_ColumnsIn2DShape(1,1) = BTTB_ColumnsIn2DShape(1,1) + alpha;
%% Creating the extended real symmetric BCCB from just the first column of the real symmetric BTTB
BCCB_Column = [];                                 % Making each block circulant, extension at each block level
for k= 1:N+1
    bttb_col =  BTTB_ColumnsIn2DShape(:,k);
    flippedColumn = flipud(bttb_col);
%     flippedColumn = [0; flippedColumn(1:N-1)];            %  for full doubled circulant even extension
    BCCB_Column = [BCCB_Column; bttb_col; flippedColumn(2:N) ];  % remove (2:N) for full circulant 
end

flippedIndexes = fliplr(2:N);                   % Making full matrix block circulant, extension at each block level
for k= 1:N-1
    bttb_col = BTTB_ColumnsIn2DShape(:,flippedIndexes(k));
    flippedColumn = flipud(bttb_col);
    BCCB_Column = [BCCB_Column; bttb_col; flippedColumn(2:N) ];
end
BCCB_ColumnsIn2DShape = reshape(BCCB_Column, [2*N, 2*N]);

IndexesBCCB = (1:2*N)';
Circulant_Indexes = [];
for i = 1:2*N
    Circulant_Indexes = [Circulant_Indexes, IndexesBCCB];
    IndexesBCCB = circshift(IndexesBCCB, 1);
end

BCCB_Matrix = []; 
for l = 1:2*N
    ColumnIndexes = Circulant_Indexes(:,l);
    CirculantBlocks = [];
    for k = 1:2*N
        CurrentBlockIndex = ColumnIndexes(k);                                           % This is the block Index
        CirculantBlock = toeplitz(BCCB_ColumnsIn2DShape(:,CurrentBlockIndex));           % This is a circulant block now
        CirculantBlocks = [ CirculantBlocks ; CirculantBlock];                           % Append Column-wise C = [C1 ; C2; C3; C2]
    end
    BCCB_Matrix  = [BCCB_Matrix, CirculantBlocks];                                        % Append Row-wise
end
%% Plotting the BTTB matrix and its eigen values
str1 = num2str(N);

figure, imagesc(BTTB_Matrix(1:N+1,1:N+1))
str2 = strcat('First Toeplitz block of the BTTB matrix,N=',str1, ', (N+1) \times (N+1)' ); 
title(str2)
colorbar

figure, imagesc(BTTB_Matrix(N+2:2*(N+1),1:N+1))
str3 = strcat('Second Toeplitz block of the BTTB matrix,N=',str1, ', (N+1) \times (N+1)' ); 
title(str3)
colorbar

figure, 
% surf(BTTB_Matrix)
imagesc( BTTB_Matrix)
str4 = strcat('BTTB Matrix, N =',str1, ', (N+1)^2 \times (N+1)^2' ); 
title(str4)
% xlabel('x')
% ylabel('y')
% zlabel('z')
%  colormap cool
colorbar
 
figure,
semilogy(flipud(eig(real(BTTB_Matrix))), 'lineWidth', 2.3)
str11 = 'Eigen values of the real symmetric BTTB matrix, ';
str2 = strcat( 'Condition Number = ', num2str(cond(BTTB_Matrix),4));
title({str11,str2})                % Using cell array
grid on
xlabel('Eigen value index')
ylabel('Magnitude of Eigen value')

%% Checking the clustering of EigenValues using the diagonal preconditioner

% NumericMatrix = BTTB_Matrix;
% iterations = 500;
% RandomVector = zeros(1, iterations); 
% Condition_Nos = zeros(1, iterations);
% 
% for k = 1: iterations
%    RandomVector(k) = rand(1,1);
%    Condition_Nos(k)= cond(NumericMatrix^(-(1-RandomVector(k)))* NumericMatrix);
% end
%  
% 
% minCond = min(Condition_Nos);
% indexes = 1:iterations;
% min_k = indexes(Condition_Nos == minCond);
% 
%     
% figure, 
% semilogy(RandomVector, Condition_Nos, '+')
% xlabel('Value of random number')
% ylabel('Condition number , \alpha random number, cond(A^{-(1-\alpha)}* A)')
% 
% minCond
% min_k 
% 
% 
% % W         = sqrt(2.5*abs(-N/2:N/2))/(N+1);             % Checking how well the fourier based preconditioner works
% % W(N/2+1)  = sqrt(1/1)/(N+1);
% % W = ones(M,1)*W ;
% % W = W.^2;
% % 
% % BTTB_ColumnsIn2DShape1 = Adjoint2DPolarDFT( W .* Compute2DPolarDFT( ImageFinal,  M ));     % Applied Diagonal operator as preconditioner
% % BTTB_ColumnsIn2DShape1 = real(BTTB_ColumnsIn2DShape1);
% % Indexes = 1:(N+1);
% % Toeplitz_Indexes = toeplitz(Indexes);   % Block Level
% % 
% % BTTB_Matrix1 = []; 
% % for l = 1:(N+1)
% %     ColumnIndexes = Toeplitz_Indexes(:,l);
% %     ToeplitzBlocks = [];
% %     for k = 1:(N+1)
% %         CurrentBlockIndex = ColumnIndexes(k);                      % This is the block Index
% %         ToeplitzBlock = toeplitz(BTTB_ColumnsIn2DShape1(:,CurrentBlockIndex));
% %         ToeplitzBlocks = [ ToeplitzBlocks ; ToeplitzBlock];    % Append Column-wise T = [T1 ; T2; T3]
% %     end
% %     BTTB_Matrix1 = [BTTB_Matrix1, ToeplitzBlocks];              % Append Row-wise
% % end
% % EigWithPreconditioner = eig(BTTB_Matrix1);
% 
% EigOriginal =  eig(NumericMatrix);
% EigWithPreconditioner = eig(NumericMatrix^(-(1-RandomVector(min_k)))* NumericMatrix);
% 
% figure, semilogx( EigOriginal , ones(1,length(EigOriginal)) , '+', real(EigWithPreconditioner) , 2* ones(1, length(EigWithPreconditioner)), '.', 'LineWidth', 2.3)
% axis([ 10^-15 10^6  0  4])
% xlabel('Eigen values')
% ylabel('Application of preconditioner')
% legend('No preconditioner',' With Preconditioner')
% title('Clustering of Eigenvalues')



%% Plotting the BCCB matrix and its eigen values comparision with FFT based and Direct approach
figure, imagesc(BCCB_Matrix)
str5 = strcat('BCCB Matrix, N =',str1, ', (2N)^2 \times (2N)^2' ); 
title(str5)
colorbar

disp('The Determinant')
det(BCCB_Matrix)

Eig_BCCB = flipud(eig(real(BCCB_Matrix)));                                              % Direct Computation of Eigen Values of a BCCB matrix
figure, 
semilogy((Eig_BCCB), 'lineWidth', 2.3)
str1 = 'Eigen/singular values of the real symmetric BCCB matrix, ';
str2 = strcat('Condition Number = ', num2str(cond(BCCB_Matrix),4));
title({str1, str2})
grid on
xlabel('Singular value index')
ylabel('Magnitude of Singular value')

FFTBasedEigenValues = reshape( fft2(BCCB_ColumnsIn2DShape) , [(2*N)^2,1]);              % FFT based eigenvalues compuations of a BCCB matrix
max(flipud(sort(real(FFTBasedEigenValues))) - Eig_BCCB)                                 % Comparison, This is exactly zero except for numerical errors in order of  10^-11 


%% Computing the inverse
figure, imagesc(inv(real(BCCB_Matrix))) 
title('Direct Inverse of the BCCB matrix ')
xlabel('x')
ylabel('y') 
colorbar
