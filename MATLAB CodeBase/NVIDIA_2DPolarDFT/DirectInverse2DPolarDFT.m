function [ ImageFinal ] = DirectInverse2DPolarDFT( Polar_Grid )
%INVERSE2DPOLARDFT is the inverse polar grid computation, given Polar_grid
%of size M x (N+1), where M is the number of angles (even)
% N+1 is the number of elements in each anglular slice of the 2D polar grid
% Solution is based on Fourier domain interpolation on circles followed by
% FrFT in radial and angular directions
% This function has been deprecated on 11/29/2015
[M, sizeN_plus1 ] = size(Polar_Grid);
N = sizeN_plus1 -1;
ImageFinal = zeros(N+1,N+1);

%% Highly Specialized Direct Solution
ImageFinal(1,1) = 1;            % Making it an impulse or a point spread function, 
% Special unique column of the BTTB matrix which is formed by product (A^H * A). This can be computed only once at the beginning for different inputs
BTTB_ColumnsIn2DShape = Adjoint2DPolarDFT( Compute2DPolarDFT( ImageFinal,  M ));     % It is (N+1) x (N+1) size output, it is the impluse response of the linear system  (A^H A)
% figure, imagesc(real(BTTB_ColumnsIn2DShape))              % Debug
% title('Impulse response or point Spread function')

BCCB_ColumnIn2DShape = makeBCCBFromBTTB_extension(BTTB_ColumnsIn2DShape,N);
% figure, imagesc(real(BCCB_ColumnIn2DShape))                % Debug
% title('Block Circulant (BCCB) extension')

b = Adjoint2DPolarDFT( Polar_Grid);                         % (N+1) x (N+1) output,  A^H y as RHS

zeroPadded_b =  padarray(b, [N-1 N-1], 'post');             % zero padded b at unknown location of reflection coefficients  extension of the input


Eigen_values_BCCB = fft2(BCCB_ColumnIn2DShape); 

c = real(ifft2(zeroPadded_b./ Eigen_values_BCCB));             % This is the right hand side of the system  " d * hat_b   = -c*padded_b "  which is "FinalImage = (A^H*A)^{-1} reflection_b"


allindexes = 1:(2*N)^2;
mysize = [2*N,2*N];
twoDIndexes = reshape(allindexes, mysize);
desiredIndexes = unique([twoDIndexes(N+2:2*N,:) , twoDIndexes(:,N+2:2*N)'], 'sorted');
c = c(desiredIndexes);          % This is the full RHS


BCCB_Matrix = createBCCBMatrixFrom2DStructure(BTTB_ColumnsIn2DShape,N) ;   % Temporary
InverseBCCBMatrix =  inv(BCCB_Matrix);

desiredRows = InverseBCCBMatrix(desiredIndexes,:);
d = desiredRows(:, desiredIndexes);    % desired columns
d= real(d);
figure, imagesc(d)
colorbar
% R = chol(d);                             %Systems like  Ax = b replaced by R?Rx = b gives x = R\(R'\b)  Assuming symmetric positive definiteness
[L,D] = ldl(d);
% hat_b = R\(R'\(-c));              %these are the reflection coefficients         
hat_b = -inv(d)*c;

reflection_b = zeroPadded_b;         % Gathering the reflection coefficients
reflection_b(desiredIndexes)  =  hat_b;

% Final inversion 
ImageFullFinal = real(ifft2(reflection_b)./ Eigen_values_BCCB);        %Solving x = (A^H*A)^-1 (A^H y) with BCCB extension   the full linear system like x = A^-1 b
ImageFinal  = ImageFullFinal(1:N+1, 1:N+1);                             % Discarding other zero pad values 



% %% Solving x = (A^H*A)^-1 A^H y   the full linear system like x = A^-1 b
% b = Adjoint2DPolarDFT( Polar_Grid);   % (N+1) x (N+1)   A^H y
% padded_b =  padarray(b, [N-1 N-1], 'post');                 % hat_b  extension of the input
% 
% Eigen_values_BCCB = fft2(BCCB_ColumnIn2DShape);  
% 
% figure, imagesc(real(Eigen_values_BCCB)/ mean(real(Eigen_values_BCCB)))
% colorbar
% 
% % % TSVD , K 
% % Phi = zeros(size(Eigen_values_BCCB));
% % tol = 5000;
% % Phi  = real(Eigen_values_BCCB) >= tol ;             % Logical array that gives only values that are greater than tolerance
% 
% % Tikhonov Regularization
% alpha = 20;                                             % Scalar Regularization
% Phi = abs(Eigen_values_BCCB).^2 / (abs(Eigen_values_BCCB).^2 + alpha^2);
% 
% SFilt = Phi ./ Eigen_values_BCCB;
% figure, imagesc(real(SFilt)/ mean(real(SFilt)))
% colorbar  
% 
% ImageFullFinal = ifft2( fft2(padded_b) .*  SFilt);
% figure, imagesc( ( real(ImageFullFinal)))
% ImageFinal  = real(ImageFullFinal(1:N+1, 1:N+1));   % We use this real() since due to overflow error small imaginary numbers may surface
 
end  

function AppendedTwoDStructure = makeBCCBFromBTTB_extension (TwoDStructure,N)               % This is just one column extension
% This function assumes that the BTTB is symmetric real matrix
% Transforming the first column of BTTB into BCCB column, requires  reflexive extentions per block as well as for entire matrix
BCCB_Column = [];                                 % Making a full block circulant with circulant blocks
for k= 1:N+1 
    bttb_col =  TwoDStructure(:,k);
    flippedColumn = flipud(bttb_col);
    BCCB_Column = [BCCB_Column; bttb_col; flippedColumn(2:N) ];            % periodic extension
end

flippedIndexes = fliplr(2:N);                   % Making full matrix block circulant
for k= 1:N-1
    bttb_col = TwoDStructure(:,flippedIndexes(k));
    flippedColumn = flipud(bttb_col);
    BCCB_Column = [BCCB_Column; bttb_col; flippedColumn(2:N) ];         % periodic extension
end

AppendedTwoDStructure = reshape(BCCB_Column, [2*N 2*N ]);
end


function BTTB_Matrix = createBTTBMatrixFrom2DStructure(BTTB_ColumnsIn2DShape,N)                % Test function useful for debugging
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
end

function BCCB_Matrix = createBCCBMatrixFrom2DStructure(BTTB_ColumnsIn2DShape,N)                  % Test function useful for debugging
BCCB_ColumnsIn2DShape = makeBCCBFromBTTB_extension (BTTB_ColumnsIn2DShape,N); 

IndexesBCCB = (1:2*N)';
Circulant_Indexes = [];
for i = 1:2*N
    Circulant_Indexes = [Circulant_Indexes, IndexesBCCB];
    IndexesBCCB = circshift(IndexesBCCB, 1);
end

count = 0;
BCCB_Matrix = []; 
for l = 1:2*N
    ColumnIndexes = Circulant_Indexes(:,l);
    CirculantBlocks = [];
    for k = 1:2*N
        CurrentBlockIndex = ColumnIndexes(k);                                           % This is the block Index
        CirculantBlock = toeplitz(BCCB_ColumnsIn2DShape(:,CurrentBlockIndex));           % This is a circulant block now
        CirculantBlocks = [ CirculantBlocks ; CirculantBlock];                           % Append Column-wise C = [C1 ; C2; C3; C2]
        
        count = count+1;
        disp(count)
    end
    BCCB_Matrix  = [BCCB_Matrix, CirculantBlocks];                                       % Append Row-wise
end
end


%% To test BCCBFrom BTTB
% 
% syms a b c d e f g h i 
% TwoDStructure = [a b c; d e f; g h i].'
% TwoDStructure =            % This is a BTTB first column in 2D shape
% [ a, d, g]
% [ b, e, h]
% [ c, f, i]
% 
% N= 2
% % Run function
% makeBCCBFromBTTB_extension (TwoDStructure,N)
% AppendedTwoDStructure =    % This is a BCCB first column in 2D shape
%  
% [ a, d, g, d]
% [ b, e, h, e]
% [ c, f, i, f]
% [ b, e, h, e]
