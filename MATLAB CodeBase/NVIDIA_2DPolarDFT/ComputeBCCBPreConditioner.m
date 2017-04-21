function [ column1In2D_BCCB ] = ComputeBCCBPreConditioner( column1In2D_BTTB )
% This function takes in first column of the real symmetric BTTB matrix
% and returns full BCCB preconditioner matrix based on T. Chan's
% preconditioning formula
% Input 1st Column of BTTB matrix which is real symmetric with symmetric blocks 
% Output the desired BCCB preconditoner 
% deprecated
[K,~] = size(column1In2D_BTTB);  % Assuming this is (N+1)^2 length

column1In2D_BCCB = zeros(size(column1In2D_BTTB));

for k = 1:K     % BC preconditioning
    column1In2D_BCCB(:,k) = Compute_TChan_Preconditioner_BC (column1In2D_BTTB(:,k));      
end

column1In2D_BCCB  = Compute_TChan_Preconditioner_BCCB (column1In2D_BCCB);  % BCCB preconditioning
% One can prepare the Complete BCCB preconditioner from the first column
end


function [BCCB_Column_2D] = Compute_TChan_Preconditioner_BCCB (columns_Circulant_2D) % BCCB preconditioning
BCCB_Column_2D = zeros(size(columns_Circulant_2D));
[K,~] = size(columns_Circulant_2D);
for k = 0:K-1
    linearCombination =  (K-k)* columns_Circulant_2D(:,k+1)+ k * columns_Circulant_2D(:,K-k );     % T. Chan's formula
    BCCB_Column_2D (:,k+1) = 1/K*linearCombination;
end
end


function [Column_Circulant] = Compute_TChan_Preconditioner_BC (column_Toeplitz) % BC preconditioning
K = length(column_Toeplitz);
Column_Circulant = zeros(K,1);
for k = 0:K-1
    linearCombination =  (K-k)* column_Toeplitz(k+1)+ k* column_Toeplitz(K-k);     % T. Chan's formula
    Column_Circulant (k+1) = 1/K*linearCombination;
end
end


% function Full_Circulant_Matrix = Make_Circulant_Matrix_FromColumn1( Column11D_Circulant)    % We dont need to make this full circulant matrix it is usually too large, but this code can be used for debugging purposes
% [K] = length(Column11D_Circulant);
% Full_Circulant_Matrix = [];
% for i = 1:K
%     Full_Circulant_Matrix = [Full_Circulant_Matrix, Column11D_Circulant ];
%     Full_Circulant_Matrix = circshift(Full_Circulant_Matrix, 1);
% end
% end