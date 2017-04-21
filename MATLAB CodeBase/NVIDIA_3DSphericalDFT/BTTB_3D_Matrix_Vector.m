function b_iterate_1D = BTTB_3D_Matrix_Vector(x_iterate_1D,BCCB_3D_EigenValues,N )
x_iterate_3D = reshape(x_iterate_1D, [N+1, N+1, N+1]);
b_iterate_3D = zeros(N+1,N+1,N+1);
toeplitzIndexes = toeplitz(1:N+1);
for l = 1:N+1
    for m = 1:N+1
        k = toeplitzIndexes(l,m);          % Toeplitz based indexes
        x_iterate_Inner_1D = reshape( x_iterate_3D(:,:,m),[(N+1)^2, 1]);
        BCCB_EigenValues = BCCB_3D_EigenValues(:,:,k);
        b_iterate_Inner_1D = BTTB_2D_Matrix_Vector(x_iterate_Inner_1D,BCCB_EigenValues,N );
        current2DImage = reshape(b_iterate_Inner_1D , [N+1,N+1]); 
        b_iterate_3D(:,:,l) = b_iterate_3D(:,:,l) + current2DImage;
    end
end
b_iterate_1D = reshape(b_iterate_3D, [(N+1)^3 1]);
end