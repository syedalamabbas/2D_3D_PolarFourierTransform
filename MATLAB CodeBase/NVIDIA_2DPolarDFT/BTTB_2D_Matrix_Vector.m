function b_iterate_1D = BTTB_2D_Matrix_Vector(x_iterate_1D,BCCB_EigenValues,N )
    inputImage       = reshape ( x_iterate_1D, [N+1, N+1]);           % Column stacking unfolding
    inputImage       = padarray(inputImage, [N-1,N-1], 'post');
    b_iterate2D      = ifft2( BCCB_EigenValues .* fft2(inputImage) ); % Matrix computations with the BCCB matrix multiplication
    b_iterate2D      = real( b_iterate2D );                           % Due to round off error the computed solutions may contain complex parts
    b_iterateRefined = b_iterate2D( 1:N+1,1:N+1 );                    % Taking just the first part of the BCCB matrix multiplication leaving last rows and cols
    b_iterate_1D     = reshape( b_iterateRefined, [(N+1)^2, 1] );     % Column stacking
end
