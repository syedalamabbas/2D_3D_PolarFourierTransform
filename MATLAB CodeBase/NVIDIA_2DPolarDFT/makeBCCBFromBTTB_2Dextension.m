function BCCB_ColumnIn2DShape = makeBCCBFromBTTB_2Dextension(BTTB_ColumnsIn2DShape,N)           % N is even and same definition as used in radon transform
% This function assumes that the BTTB is symmetric real matrix
BCCB_ColumnIn2DShape = zeros(2*N, 2*N);                                 % Making a full block circulant with circulant blocks
for k= 1:N+1
    bttb_col =  BTTB_ColumnsIn2DShape(:,k);
    flippedColumn = flipud(bttb_col);
    BCCB_ColumnIn2DShape(:,k) = [bttb_col; flippedColumn(2:N)];         % periodic extension
end
flippedIndexes = fliplr(2:N);                    % Making full matrix block circulant
for k= 1:N-1
    bttb_col = BTTB_ColumnsIn2DShape(:,flippedIndexes(k));
    flippedColumn = flipud(bttb_col);
    BCCB_ColumnIn2DShape(:,k+N+1) = [bttb_col; flippedColumn(2:N)];     % periodic extension
end
end