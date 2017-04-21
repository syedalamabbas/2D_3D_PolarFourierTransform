clear ;
clc;

Size = 7;

%% Symbols for my Image
Image = sym('a',[Size Size]);


%% Symbols for scaling
N = Size -1;
M = (Size+1);                                    % Number of angles
L = (M-2)/4;                                   % Number of levels
if(rem(M-2,4) ~= 0)
    L = ceil (L);                % use + 1 to compute for 45
end

PolarSymbolic = sym('p', [M, Size]);
assume(PolarSymbolic,'real');

alphas = sym('alpha_',[1 L]);
assume(alphas,'real');
betas = sym('beta_',[1 L]);
assume(betas ,'real');
syms one
assume(one ,'real');


%%  X-Butterfly & Y-Butterfly
% syms FrFT_XAxis_Uniform FrFT_YAxis_Uniform  PolarGrid PolarGridStackedColumnWise
% for p = 1: Size      % For row
%     for q = 1: Size   % For columns
%         FrFT_XAxis_Uniform(p,q) = SymbolicFrFT_Uniform( Image(p,:), alphas(1), q-1-(Size-1)/2, Size );
%         %         simplify (FrFT_XAxis_Uniform(r,c))
%     end
% end
%
% %% computing for beta
% spacing = [-(Size-1)/2 : (Size-1)/2];
% IndexBetas = spacing' * spacing ;
% Line_1 = sum ( FrFT_XAxis_Uniform .* exp(-2*1i*pi*IndexBetas* betas(1)/ Size) )
%
%
% % Line_1(1)
% % Line_1(2)
% % Line_1(3)
% % Line_1 (6)
%
%
% %% Preparing the matrix for ones and zeros

% %% Solving for first line first element testing
% % Get First Row for First element
% % A x = b , where b = Line_1(1), x is column stacked image
%
% for h = 1:1 % Size
% for z = 1: Size*Size
%     A(1,z) = subs(Line_1(h), columnStackedImage (indexes (indexes ~= z)), zerosCell ) ;
% end
%
% InnerRow1 = (subs(A(1,:),columnStackedImage, onesCell));
% InnerRow1 = (simplify(InnerRow1));
% % InnerRow1
% logOfInnerRow1 = (simplify(log(InnerRow1),'IgnoreAnalyticConstraints', true));
% factoredLogInnerRow = logOfInnerRow1*Size/(4*1i*pi);
% factoredLogInnerRow
%
% % columnStackedImage
%
% end

%% Now doing it for the entire polar grid
columnStackedImage = reshape(Image,[1,Size*Size]);
[ PolarGrid ] = ComputeSymbolic2DPolarDFT( Image,  M , alphas, betas, one);

% rowStackedPolarVals = reshape ( PolarSymbolic .', [1, M*Size]);                                        % M(N+1) x 1
% [ ImageAdjoint ] = ComputeSymbolicAdjoint2DPolarDFT( PolarSymbolic, alphas, betas, one_alpha );
%
% %% Evaluate the fully numeric PolarGrid to verify symbolic operations are correct, check all entries manually by evaluating each row or column of Polar Grid
%  alpha_Values  = num2cell(cosd((1:L)* 180/M));
%  beta_Values = num2cell(sind((1:L)* 180/M));

% TestImage = rand(Size);
% myImage = reshape(TestImage, [1 Size*Size]);
%
% PolarGrid = subs(PolarGrid, alphas, alpha_Values);
% PolarGrid = subs(PolarGrid, betas, beta_Values);
% PolarGrid = subs (PolarGrid, one_alpha, 1);
% PolarGrid = subs(PolarGrid, columnStackedImage, num2cell(myImage));
%
% DirectNumericSolutionForward = Compute2DPolarDFT( TestImage, M );  % the symbolic computations should give us the same values when substituted
%
% for k = 1: 4
%     index = randi(M);
%     PolarGrid(index,1:2)
%     DirectNumericSolutionForward (index, 1:5)
% end

%% Testing the numeric and symbolic adjoints
% TestPolarGrid = rand(M, Size) + 1i* rand(M,Size);  % generating random complex numbers
% myPolarSolution = reshape (TestPolarGrid.', [1 M*Size]);
%
% DirectNumericSolutionReverse = Adjoint2DPolarDFT(TestPolarGrid);
% ImageAdjoint = subs(ImageAdjoint, alphas, alpha_Values);
% ImageAdjoint = subs(ImageAdjoint, betas, beta_Values);
% ImageAdjoint = subs(ImageAdjoint, one_alpha, 1);
% ImageAdjoint = subs(ImageAdjoint, rowStackedPolarVals, num2cell(myPolarSolution) );
%
% for k = 1: 4
%     index = randi(Size);
%     ImageAdjoint(index,1:2)
%     DirectNumericSolutionReverse (index, 1:5)
% end


%% iterating through all lines and forging N^2 x N^2 matrix Ax = y where y is the stacked polar lines from 0 to 180

zerosCell = cell(1, Size*Size-1 );
for h = 1: Size*Size-1
    zerosCell{h} = 0;
end
onesCell = cell(1,Size*Size );
for h = 1: Size*Size
    onesCell{h} = 1;
end
indexes = 1:Size*Size;

syms A Temp  % the final A how it should look like

counter = 1;
for m = 1: M  % Iterate Through all angles
    CurrentLine  = PolarGrid(m,:);
    %         CurrentLine
    for h = 1: Size % Iterate through all elements of a line
        
        parfor z = 1: Size^2                                 % Special SpeedUP ???
            Temp(z) = subs(CurrentLine(h), columnStackedImage (indexes (indexes ~= z)), zerosCell ) ;
        end
        
        InnerRow1 = (subs(Temp,columnStackedImage, onesCell));
        %                 InnerRow1 = (simplify(InnerRow1));
        %                 logOfInnerRow1 = (simplify(log(InnerRow1),'IgnoreAnalyticConstraints', true));
        %                 factoredLogInnerRow = logOfInnerRow1*Size/(1i*pi);
        fprintf('Just finished iteration #%d\n', counter);
        
        A (counter,1: Size^2) =InnerRow1; %factoredLogInnerRow ;    % Skip dc value
        counter = counter+1;
    end
end

disp('Simplifying just A')
tic
A = simplify(A);
toc
disp('Simplified just A!')

hermA = A';
disp('Computed the Hermitian A(A^H)')


% zerosCell = cell(1, M*Size-1 );
% for h = 1: M*Size-1
%     zerosCell{h} = 0;
% end
% onesCell = cell(1,M*Size );
% for h = 1: M*Size
%     onesCell{h} = 1;
% end
% indexes = 1:M*Size;
%
% syms A_H Temp  % the final A how it should look like
%
% counter = 1;
% for m = 1: Size  % Iterate Through all angles
%     %      if (m ~= zeroIndex )     % Skip 0 degrees and 90 lines
%     CurrentLine  = ImageAdjoint(m,:);
%     %         CurrentLine
%     for h = 1: Size % Iterate through all elements of a line
%         %             if (h ~= (Size-1)/2+1)
%         for z = 1: Size*M
%             Temp(z) = subs(CurrentLine(h), rowStackedPolarVals (indexes (indexes ~= z)), zerosCell ) ;
%         end
%
%         InnerRow1 = (subs(Temp,rowStackedPolarVals, onesCell));
%         %                 InnerRow1 = (simplify(InnerRow1));
%         %                 logOfInnerRow1 = (simplify(log(InnerRow1),'IgnoreAnalyticConstraints', true));
%         %                 factoredLogInnerRow = logOfInnerRow1*Size/(1i*pi);
%         fprintf('Just finished iteration #%d\n', counter);
%
%         A_H (counter,1: Size*M) =InnerRow1; %factoredLogInnerRow ;    % Skip dc value
%         counter = counter+1;
%         %          end
%     end
%     %        end
% end
% A_H = simplify(A_H)
%
% hermA - A_H % This should be zero for the Adjoint to be correct


product = hermA* A  ;   % A^H * A          % permutations of the polar grid w.r.t row or anglular slices doesn't affect the end final product
disp('Computed the product')

% memory 
% disp('Simplifying the product, this may take time')
% tic
% for m = 1: (N+1)^2
%     parfor k = 1: (N+1)^2
%         product(k,m) = simplify(product(k,m), 'IgnoreAnalyticConstraints', true);
%         disp(['Simplification Iteration #' , num2str( (m-1)*(N+1)^2 + k )] )
%     end
% end
% toc
% disp('Simplified the full product!')

column1  = product(:,1);
disp('Simplifying the first Column only, this may take time')
tic
for k = 1: (N+1)^2
    column1(k) = real(simplify(column1(k), 'IgnoreAnalyticConstraints', true ));
    disp(['Simplification Iteration #' , num2str(  k )] )
end
toc
disp('Simplified the first column of real symmetric block toeplitz with toeplitz blocks!')


% prodReal = real(product);
% prodImag = imag(product);
% disp('Computed the split  product')



T = sym('T_dirc_', [1 N+1]);
counter = 0;
disp('These are the individual Toeplitz Blocks Matrix columns !')
for k = 1:(N+1)
    offset = counter*(N+1) ;
    disp(T(k) )
    disp(column1( (1+ offset):((N+1)+offset)))
    counter = counter +1;
end

%% Manually enlisting for special case
% if((N+1)== 11)
%     syms T_dirc_1 T_dirc_2 T_dirc_3 T_dirc_4 T_dirc_5 T_dirc_6 T_dirc_7 T_dirc_8 T_dirc_9 T_dirc_10 T_dirc_11
%     counter = 0;
%     offset = counter*(N+1)
%     T_dirc_1 = ( prodReal(1:(N+1), (1+ offset):((N+1)+offset)));
%     counter = counter +1;
%     offset = counter*(N+1)
%     T_dirc_2 = ( prodReal(1:(N+1), (1+ offset):((N+1)+offset)));
%     counter = counter +1;
%     offset = counter*(N+1)
%     T_dirc_3 = ( prodReal(1:(N+1), (1+ offset):((N+1)+offset)));
%     counter = counter +1;
%     offset = counter*(N+1)
%     T_dirc_4 = ( prodReal(1:(N+1), (1+ offset):((N+1)+offset))); 
%     counter = counter +1;
%     offset = counter*(N+1)
%     T_dirc_5 = ( prodReal(1:(N+1), (1+ offset):((N+1)+offset)));
%     counter = counter +1;
%     offset = counter*(N+1)
%     T_dirc_6 = ( prodReal(1:(N+1), (1+ offset):((N+1)+offset)));
%     counter = counter +1;
%     offset = counter*(N+1)
%     T_dirc_7 = ( prodReal(1:(N+1), (1+ offset):((N+1)+offset)));
%     counter = counter +1;
%     offset = counter*(N+1)
%     T_dirc_8 = ( prodReal(1:(N+1), (1+ offset):((N+1)+offset)));
%     counter = counter +1;
%     offset = counter*(N+1)
%     T_dirc_9 = ( prodReal(1:(N+1), (1+ offset):((N+1)+offset)));
%     counter = counter +1;
%     offset = counter*(N+1)
%     T_dirc_10 = ( prodReal(1:(N+1), (1+ offset):((N+1)+offset)));
%     counter = counter +1;
%     offset = counter*(N+1)
%     T_dirc_11 = ( prodReal(1:(N+1), (1+ offset):((N+1)+offset)));
%     
% end

%% Numeric Matrix product
prodReal = real(product);
alpha_Values  = (cosd((1:L)* 180/M));
beta_Values = (sind((1:L)* 180/M));
NumericProdReal = subs(prodReal, alphas, alpha_Values );
NumericProdReal = subs(NumericProdReal,  betas, beta_Values );
NumericProdReal = subs(NumericProdReal, one, 1 );
NumericProdReal = eval(NumericProdReal)


NumericA = subs(A, alphas, alpha_Values );
NumericA = subs(NumericA,  betas, beta_Values );
NumericA = subs(NumericA, one, 1 );
NumericA = eval(NumericA)



%% Verify if the column can be retrieved without explicitly constructing the matrix A^H A
ImpulseInput = zeros(N+1,N+1); 
ImpulseInput(1,1) = 1;
ImpulseResponse = Adjoint2DPolarDFT( Compute2DPolarDFT( ImpulseInput,  M ));
% LHS = reshape( NumericProdReal(:,1), [N+1, N+1])
% RHS = real(ImpulseResponse)                           % Some small imaginary components maybe present, ignore 
% LHS - RHS                                             % Should be very small of order of 10^-12
%  

%% Creating Reflection BTTB system as defined
% dash_A = sym('dummy', size(A));
% rowOffset = 1;
% newA = A;
% for k = 1: M
%     ElementsInStrip = newA ( rowOffset : rowOffset+ N,:);
% 
%     colOffset = 1;
%     for h=1:N+1
%          SubElementsInStrip2D = ElementsInStrip(:, colOffset : colOffset+N);                % (N+1) x (N+1) block
% %          CircShift1 = (circshift(SubElementsInStrip2D,2,2));                     % Shift along rows
% %          CircShift2 = (circshift(CircShift1,2,1));                               % Shift alone columns
%          ElementsInStrip(:, colOffset : colOffset+N) = SubElementsInStrip2D;
%          colOffset  = colOffset + N+1;
%     end
%     
% %      ElementsInStrip = circshift(ElementsInStrip,2,2);                     % Full strip
%     
%     dash_A(rowOffset : rowOffset+N, :) =  ElementsInStrip;
%     rowOffset = rowOffset + N+1;
% end

syms l m n o l_conj m_conj n_conj o_conj
sym_herm_A = [ 
         0 0 0 l      1 l_conj 0 0 0 m      1 m_conj; 
         0 0 0 n      1 n_conj 0 0 0 n_conj 1 n;     
         0 0 0 m_conj 1 m      0 0 0 l_conj 1 l; 
         0 0 0 o      1 o_conj 0 0 0 o      1 o_conj; 
         0 0 0 1      1 1      0 0 0 1      1 1;
         0 0 0 o_conj 1 o      0 0 0 o_conj 1 o ;
         0 0 0 m      1 m_conj 0 0 0 l      1 l_conj ; 
         0 0 0 n_conj 1 n      0 0 0 n      1 n_conj; 
         0 0 0 l_conj 1 l      0 0 0 m_conj 1 m; 
         ];
     
sym_A = sym_herm_A';

sym_ProductNew = sym_herm_A* sym_A;
     
% sym_ProductNew = subs (sym_ProductNew, conj(l),l_conj );
% sym_ProductNew = subs (sym_ProductNew, conj(l_conj),l );
% sym_ProductNew = subs (sym_ProductNew, conj(m),m_conj );
% sym_ProductNew = subs (sym_ProductNew, conj(m_conj),m );
% sym_ProductNew = subs (sym_ProductNew, conj(n),n_conj );
% sym_ProductNew = subs (sym_ProductNew, conj(n_conj),n );
% sym_ProductNew = subs (sym_ProductNew, conj(o),o_conj );
% sym_ProductNew = subs (sym_ProductNew, conj(o_conj),o );


sym_ProductNew

