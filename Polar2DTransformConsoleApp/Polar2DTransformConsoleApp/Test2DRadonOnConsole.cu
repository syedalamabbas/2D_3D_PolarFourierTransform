#include <arrayfire.h>
#include <af/util.h>
#include <cuda_runtime.h>

using namespace af;

array cexp(const array &in)
{
	if (!in.iscomplex()) return exp(in);
	return exp(real(in))*complex(cos(imag(in)), sin(imag(in)));
}

static const double elements_Image[] =
{
	0.1117, 0.5606, 0.6126, 0.8452, 0.8908, 0.4843, 0.4896, 0.5386, 0.2703, 2, 2,
	0.1363, 0.9296, 0.9900, 0.7386, 0.9823, 0.8449, 0.1925, 0.6952, 0.2085, 2, 2,
	0.6787, 0.6967, 0.5277, 0.5860, 0.7690, 0.2094, 0.1231, 0.4991, 0.5650, 2, 2,
	0.4952, 0.5828, 0.4795, 0.2467, 0.5814, 0.5523, 0.2055, 0.5358, 0.6403, 2, 2,
	0.1897, 0.8154, 0.8013, 0.6664, 0.9283, 0.6299, 0.1465, 0.4452, 0.4170, 2, 2,
	0.4950, 0.8790, 0.2278, 0.0835, 0.5801, 0.0320, 0.1891, 0.1239, 0.2060, 2, 2,
	0.1476, 0.9889, 0.4981, 0.6260, 0.0170, 0.6147, 0.0427, 0.4904, 0.9479, 2, 2,
	0.0550, 0.0005, 0.9009, 0.6609, 0.1209, 0.3624, 0.6352, 0.8530, 0.0821, 2, 2,
	0.8507, 0.8654, 0.5747, 0.7298, 0.8627, 0.0495, 0.2819, 0.8739, 0.1057, 2, 2,
	0.5386, 0.2703, 0.1420, 0.8604, 0.0309, 0.5590, 0.1182, 0.8182, 0.9052, 2, 2,
	0.6952, 0.2085, 0.1665, 0.9344, 0.9391, 0.8541, 0.9884, 0.1002, 0.6754, 2, 2
}; // Additional lines for 45 degrees test

void SplitMultiplyComplex(array& A, array&  B, array& realRealPart, array& realImagPart, array& imagRealPart, array& imagImagPart)
{
	// Consider multiplication of complex numbers A = (a+ib); B = (c+id)	
	//array A, B;
	//array realRealPart;             // ac
	//array realImagPart;				// ad
	//array imagRealPart;				// bc
	//array imagImagPart;				// bd
	realRealPart = real(A)*real(B);
	realImagPart = real(A)*imag(B);
	imagRealPart = imag(A)*real(B);
	imagImagPart = imag(A)*imag(B);
}

array Compute2DColumnwise_FrFTUniform(array & Image2D, array& ColumnScales_1D, int d_NoOfElements, int d_NoOfScales)
{
	int debugFlag_ToPrint = 0;
	
	/*-----------------------------------Preparing Padded & Tiled Imag2D --------------------------------------------*/
	array Zeros = constant(0, d_NoOfElements, d_NoOfElements, f32);
	array Zero_Padded_Image2D = join(0, Image2D, Zeros);
	//af_print(Zero_Padded_Image2D);
	array Image2D_Tiled = tile(Zero_Padded_Image2D, 1, 1, d_NoOfScales);
	//af_print(Image2D_Tiled);

	//printf("Successfully prepared the padded and tiled images of size %d-by-%d-by-%d  !", 2 * d_NoOfElements, d_NoOfElements, d_NoOfScales);
	//printf("\n");

	/*-------------------------------------------Creating Index Cubes and Sequences----------------------------------------------------*/
	int N = d_NoOfElements - 1;
	array leftSideIndexes = array(seq(0, N, 1));
	array rightSideIndexesOnes = -1 * array(seq(1, d_NoOfElements, 1));
	array rightSideIndexesZeros = constant(0, d_NoOfElements, 1);
	array rightSideIndexesN_2 = constant(N / 2, d_NoOfElements, 1);

	if (debugFlag_ToPrint)
	{
		af_print(leftSideIndexes);
		//af_print(leftSideIndexes.isdouble());
		//af_print(flip( rightSideIndexesOnes,0));
		//af_print(rightSideIndexesZeros);
		//af_print(rightSideIndexesN_2);
	}


	array indexedElementsEn = join(0, leftSideIndexes, flip(rightSideIndexesOnes, 0));
	array indexedElementsPre = join(0, leftSideIndexes, rightSideIndexesN_2);     /* This is for Keeping pre and post multiplication factor upper half only*/
	array indexedElementsPost = join(0, leftSideIndexes, rightSideIndexesZeros);
	if (debugFlag_ToPrint)
	{
		//af_print(indexedElementsEn);
		//af_print(indexedElementsPre_Post);
	}


	array indexedElements_Tiled_En = tile(pow(indexedElementsEn, 2), 1, d_NoOfElements, d_NoOfScales);
	array indexedElements_Tiled_PreMulti = tile(indexedElementsPre - N / 2, 1, d_NoOfElements, d_NoOfScales);
	array indexedElements_Tiled_PostMulti = tile(indexedElementsPost, 1, d_NoOfElements, d_NoOfScales);
	
	if (debugFlag_ToPrint)
	{
		af_print(indexedElements_Tiled_En.col(0));
		af_print(indexedElements_Tiled_PreMulti.col(0));
		af_print(indexedElements_Tiled_PostMulti.col(0));
	}


	/*--------------------------Creating FrFT scale cubes------------------------------------*/
	array ColumnScales_1D_FullTiled;
	if (d_NoOfScales == 1)
	{
		//af_print(ColumnScales_1D);
		ColumnScales_1D_FullTiled = constant(ColumnScales_1D.scalar<float>(), 2 * d_NoOfElements, d_NoOfElements);
		//af_print(ColumnScales_1D_FullTiled);
	}
	else
	{
		array ColumnScales_1D_Mods  = moddims(ColumnScales_1D, 1, 1, d_NoOfScales);
		array ColumnScales_1D_Tiled_depth = tile(ColumnScales_1D_Mods, 2 * d_NoOfElements, d_NoOfElements, 1);
		ColumnScales_1D_FullTiled = moddims(ColumnScales_1D_Tiled_depth, 2 * d_NoOfElements, d_NoOfElements, d_NoOfScales);
	}
	
	if (debugFlag_ToPrint)
	{
		//af_print(ColumnScales_1D_Mods);
		//af_print(ColumnScales_1D_Tiled_depth);
		//af_print(ColumnScales_1D_FullTiled);
	}


	/*-------------------Precomputing the Essential Sequence Cubes :: All complex values here --------------------*/
	/*array imaginaryUnit_Tiled = tile(i, 2 * d_NoOfElements, d_NoOfElements, d_NoOfScales);
	if (debugFlag_ToPrint)
	{
		af_print(imaginaryUnit_Tiled.slice(0))
	}*/
	array Seq_En = cexp(complex(0, -af::Pi * indexedElements_Tiled_En * ColumnScales_1D_FullTiled / d_NoOfElements));   /* E(n) as defined in the paper*/

	//array temp = array(seq(0, 5)).as(f64);
	//af_print(temp);
	//af_print(complex(0,temp));
	/*af_print(cexp(-complex(0,af::Pi * indexedElements_Tiled_En * ColumnScales_1D_FullTiled / d_NoOfElements)));
	array Seq_En = cexp(-imaginaryUnit_Tiled*af::Pi * indexedElements_Tiled_En * ColumnScales_1D_FullTiled / d_NoOfElements);
	af_print(Seq_En);*/


	array Ones = constant(1, d_NoOfElements, d_NoOfElements, f32);
	array subtractValues = tile(join(0, Zeros, Ones), 1, 1, d_NoOfScales);			/* This is for  Keeping pre and post multiplication factor upper half only*/
	array PreMultiplicationFactor = cexp(complex (0,  af::Pi * indexedElements_Tiled_PreMulti * ColumnScales_1D_FullTiled * N / d_NoOfElements)) - subtractValues;
	array PostMultiplicationFactor = cexp(complex(0, af::Pi *  indexedElements_Tiled_PostMulti * ColumnScales_1D_FullTiled * N / d_NoOfElements)) - subtractValues;

	if (debugFlag_ToPrint)
	{
		af_print(Seq_En);
		af_print(PreMultiplicationFactor.col(0));
		af_print(PostMultiplicationFactor.col(0));
	}


	/*--------------------Preprocessing Cubes-----------------------*/
	array Image2D_Tiled_PreMulti = Image2D_Tiled * PreMultiplicationFactor;
	array Image2D_Tiled_PreMulti_SeqEn = Image2D_Tiled_PreMulti * Seq_En;

	if (debugFlag_ToPrint)
	{
		af_print(Image2D_Tiled_PreMulti);
		af_print(Seq_En);
		af_print(Image2D_Tiled_PreMulti_SeqEn);
	}


	/*-------------------Computing Convolution--------------------*/
	array firstFFT_X = fft(Image2D_Tiled_PreMulti_SeqEn);
	array secondFFT_X = fft(conjg(Seq_En));
	array interim_FrFT_X = ifft(firstFFT_X * secondFFT_X);

	if (debugFlag_ToPrint)
	{
		af_print(firstFFT_X.col(0));
		af_print(secondFFT_X.col(0));
		af_print(interim_FrFT_X.cols(0, N));
	}



	/*-------------------Postprocessing-----------------------------*/
	array  FrFT_Image_X = interim_FrFT_X * Seq_En * PostMultiplicationFactor;

	/*--------------------Grab only the top half elements drop overlapping------------------*/
	array FrFT_Image_X_Cube = FrFT_Image_X.rows(0, N);

	if (debugFlag_ToPrint)
	{
		af_print(FrFT_Image_X);
		af_print(FrFT_Image_X_Cube.cols(0, N));
	}

	return FrFT_Image_X_Cube;
}

void Compute2DComplementaryLines_FrFTVariableScales(array & OneD_FrFT, array& ColumnScales_1D, array& final2DFrFTImage, array& final2DFrFTConjImage, int d_NoOfElements, int d_NoOfScales)
{
	int debugFlag_ToPrint = 0;

	int N = d_NoOfElements - 1;
	array lineSpacing = array(seq(-N / 2, N / 2));
	array lineSpacing_tiled_Y = tile(lineSpacing, 1, d_NoOfElements);
	array lineSpacing_Square = lineSpacing_tiled_Y * lineSpacing_tiled_Y.T();

	array lineSpacing_Square_TiledLevel = tile(lineSpacing_Square, 1, 1, d_NoOfScales);
	//printf("\nCreated preliminary indexed matrices and sequences, E_n, PreMultiplication, PostMultiplication and Linespacing .\n");

	af_print(lineSpacing_Square);
	if (debugFlag_ToPrint)
	{
		//af_print(lineSpacing);
		//af_print(lineSpacing_tiled_Y);
		af_print(lineSpacing_Square);
		af_print(lineSpacing_Square_TiledLevel);
	}

	array beta_Tiled;
	if (d_NoOfScales == 1)
	{
		beta_Tiled = constant(ColumnScales_1D.scalar<float>(), d_NoOfElements, d_NoOfElements);
	}
	else
	{ 
		array beta_Mods = moddims(ColumnScales_1D, 1, 1, d_NoOfScales);
		array beta_Tiled_depth = tile(beta_Mods, d_NoOfElements, d_NoOfElements, 1);
		beta_Tiled = moddims(beta_Tiled_depth, d_NoOfElements, d_NoOfElements, d_NoOfScales);
	}
	//printf("\nCreated alpha and beta matrices.\n");

	if (debugFlag_ToPrint)
	{
		//af_print(beta_Tiled_depth);
		af_print(beta_Tiled);
	}


	/*-------------------Precomputing the Essential Sequence Cubes :: All complex values here --------------------*/
	array BetaFactor = cexp( complex(0, -2 * af::Pi * lineSpacing_Square_TiledLevel * beta_Tiled / d_NoOfElements));


	// Consider multiplication of complex numbers A = (a+ib); B = (c+id)	
	array realRealPart;             // ac
	array realImagPart;				// ad 
	array imagRealPart;				// bc
	array imagImagPart;				// bd

	SplitMultiplyComplex(OneD_FrFT, BetaFactor, realRealPart, realImagPart, imagRealPart, imagImagPart);

	array tempSeq_X = sum(complex(realRealPart - imagImagPart, realImagPart + imagRealPart)); // sum(flip(FrFT_Image_X_Cube, 0) *BetaFactor);
	
	array tempSeqConj_X = sum(complex(realRealPart + imagImagPart, imagRealPart - realImagPart));//  sum(flip(FrFT_Image_X_Cube, 0) * conjg(BetaFactor));
	
	if (debugFlag_ToPrint)
	{
		af_print(BetaFactor.cols(0, N));

		af_print(tempSeq_X);
		//af_print(sum(flip(FrFT_Image_X_Cube, 0) *BetaFactor));
		af_print(tempSeqConj_X);
		//af_print(sum(flip(FrFT_Image_X_Cube, 0) * conjg(BetaFactor)));
	}

	final2DFrFTImage = moddims(tempSeq_X, d_NoOfElements, d_NoOfScales).T();
	final2DFrFTConjImage = moddims(tempSeqConj_X, d_NoOfElements, d_NoOfScales).T();
}

void FullyVectorized2DPolarTransform(array& Image, int d_NoOfElements, int d_NoOfAngles, int d_NoOfLevels, array& FinalPolarGrid)
{
	int debugFlag_ToPrint = 0;
	int N = d_NoOfElements - 1; /* it is always even as described in the paper*/

	//array Image = randu(d_NoOfElements, d_NoOfElements, f64);
	printf("\nCreating a %d-by-%d elements of an ArrayFire object, it has %d levels.", d_NoOfElements, d_NoOfElements, d_NoOfLevels);
	printf("\n");

	if (debugFlag_ToPrint)
	{
		af_print(Image);
		af_print(Image.row(0));
		af_print(Image.col(0));
	}

	/*--------------------------Creating Alpha Cubes------------------------------------*/
	array alpha_Levels = cos(af::Pi / d_NoOfAngles * array(seq(1, d_NoOfLevels)));
	if (debugFlag_ToPrint)
	{
		af_print(alpha_Levels);
	}
	/*--------------------------Creating Beta Cubes---------------------------------------*/
	array beta_Levels = sin(Pi / d_NoOfAngles * array(seq(1, d_NoOfLevels)));
	//af_print(beta_Levels);
	//printf("\nCreated alpha and beta matrices.\n");

	/*-------------------Precomputing the Essential Sequence Cubes :: All complex values here --------------------*/
	array lineSpacing = array(seq(-N / 2, N / 2));
	array lineSpacing_tiled_Y = tile(lineSpacing, 1, d_NoOfElements);
	array lineSpacing_Square = lineSpacing_tiled_Y * lineSpacing_tiled_Y.T();
	array ZeroNinty_Factor = cexp(complex(0, -2 * af::Pi * lineSpacing_Square * 1 / d_NoOfElements));
	if (debugFlag_ToPrint)
	{
		//af_print(lineSpacing);
		//af_print(lineSpacing_tiled_Y);
		af_print(lineSpacing_Square);
		af_print(ZeroNinty_Factor);
	}



	/*-------------------- First dimension uniform FrFT for each Image per level-----------------------*/
	array FrFT_Image_X_Cube = Compute2DColumnwise_FrFTUniform(Image.T(), alpha_Levels, d_NoOfElements, d_NoOfLevels);
	array FrFT_Image_Y_Cube = Compute2DColumnwise_FrFTUniform(Image, alpha_Levels, d_NoOfElements, d_NoOfLevels);

	FrFT_Image_X_Cube = FrFT_Image_X_Cube.T();       // Now it needs operation to the other dimension
	FrFT_Image_Y_Cube = FrFT_Image_Y_Cube.T();


	/*--------------------Finally all computations for  the  Polar Grid-----------*/
	//   Computing for all the grid expect two special indexes
	array levelSeq = array(seq(0, d_NoOfLevels - 1));

	array finalIndexSeq1_X = 1 + levelSeq;

	array finalIndexSeq2_X = d_NoOfAngles - finalIndexSeq1_X;

	array finalIndexSeq3_Y = d_NoOfAngles / 2 - finalIndexSeq1_X;

	array finalIndexSeq4_Y = d_NoOfAngles / 2 + finalIndexSeq1_X;

	if (debugFlag_ToPrint)
	{
		af_print(levelSeq);
		af_print(finalIndexSeq1_X);
		af_print(finalIndexSeq2_X);
		af_print(finalIndexSeq3_Y);
		af_print(finalIndexSeq4_Y);
	}

	//array A = constant(1, 10);
	//array B = A; // B and A point to same memory
	//float *d_A = A.device<float>();
	//cudaMemset(d_A, 0, A.bytes());
	//af_print(B); //  all zeros since pointed to same memory
	//array element = finalIndexSeq3_Y.row(0);
	//af_print(element);
	//float rho = element.scalar<float>();
	//printf("%f", rho);

	array finalSeq_X, finalSeqConj_X;
	Compute2DComplementaryLines_FrFTVariableScales(FrFT_Image_X_Cube, beta_Levels, finalSeq_X, finalSeqConj_X, d_NoOfElements, d_NoOfLevels);
	finalSeqConj_X = flip(finalSeqConj_X, 1);             // Special operation
	if (debugFlag_ToPrint)
	{
		af_print(finalSeq_X);
		af_print(finalSeqConj_X);
	}

	array finalSeq_Y, finalSeqConj_Y;
	Compute2DComplementaryLines_FrFTVariableScales(FrFT_Image_Y_Cube, beta_Levels, finalSeq_Y, finalSeqConj_Y, d_NoOfElements, d_NoOfLevels);


	// Removing just 2 redundant computations for 45 degree incase
	if (0 == remainder(d_NoOfAngles, 4))
	{
		finalIndexSeq3_Y = finalIndexSeq3_Y.rows(0, d_NoOfLevels - 2);          // Removing just the last rows from 4 structures
		finalSeq_Y = finalSeq_Y.rows(0, d_NoOfLevels - 2);
		finalIndexSeq4_Y = finalIndexSeq4_Y.rows(0, d_NoOfLevels - 2);
		finalSeqConj_Y = finalSeqConj_Y.rows(0, d_NoOfLevels - 2);
	}

	if (debugFlag_ToPrint)
	{
		af_print(finalIndexSeq3_Y);
		af_print(finalIndexSeq4_Y);
		af_print(finalSeq_Y);
		af_print(finalSeqConj_Y);
	}


	//   Computing seperately for two special indexes
	float zeroIndex = 0;
	float nintyIndex = d_NoOfAngles / 2;
	float values[] = { zeroIndex, nintyIndex };
	array SpecialTwoIndexes(2, 1, values);


	array ZeroLineFrFT_Image_X_Cube = FrFT_Image_Y_Cube.slice(zeroIndex).col(N / 2);

	array NintyLineFrFT_Image_Y_Cube = FrFT_Image_X_Cube.slice(zeroIndex).col(N / 2);


	if (debugFlag_ToPrint)
	{
		af_print(FrFT_Image_Y_Cube);
		af_print(ZeroLineFrFT_Image_X_Cube);
		//af_print(FrFT_Image_X_Cube);
		af_print(NintyLineFrFT_Image_Y_Cube);
	}


	array DFTZeroLine = sum(tile(ZeroLineFrFT_Image_X_Cube, 1, d_NoOfElements) *ZeroNinty_Factor);

	array DFTNinetyLine = sum(tile(NintyLineFrFT_Image_Y_Cube, 1, d_NoOfElements) *ZeroNinty_Factor);

	array SpecialTwoLines = join(0, DFTZeroLine, DFTNinetyLine);

	if (debugFlag_ToPrint)
	{
		af_print(DFTZeroLine);
		af_print(DFTNinetyLine);
		af_print(SpecialTwoLines);
	}

	array UnsortedIndexes = join(0, join(0, join(0, join(0, finalIndexSeq1_X, finalIndexSeq2_X), finalIndexSeq3_Y), finalIndexSeq4_Y), SpecialTwoIndexes);
	array tiledUnsortedIndexes = tile(UnsortedIndexes, 1, d_NoOfElements);
	array UnsortedPolarGrid = join(0, join(0, join(0, join(0, finalSeq_X, finalSeqConj_X), finalSeq_Y), finalSeqConj_Y), SpecialTwoLines);

	if (debugFlag_ToPrint)
	{
		af_print(tiledUnsortedIndexes);
		af_print(UnsortedPolarGrid);
	}

	array FinalPolarGridReal;// = constant(0, d_NoOfElements, d_NoOfAngles, c64);
	array Output_Keys_Sorted;
	sort(Output_Keys_Sorted, FinalPolarGridReal, tiledUnsortedIndexes, real(UnsortedPolarGrid));
	if (debugFlag_ToPrint)
	{
		af_print(Output_Keys_Sorted);
		af_print(FinalPolarGridReal);
	}

	array FinalPolarGridImag;// = constant(0, d_NoOfElements, d_NoOfAngles, c64);
	array Output_Keys_Sorted2;
	sort(Output_Keys_Sorted2, FinalPolarGridImag, tiledUnsortedIndexes, imag(UnsortedPolarGrid));
	
	// Put it in the output
	FinalPolarGrid = complex(FinalPolarGridReal, FinalPolarGridImag);
	if (debugFlag_ToPrint)
	{
		af_print(FinalPolarGrid);
	}
	
}

void VectorizedPerLevel2DPolarTransform(array& Image, int d_NoOfElements, int d_NoOfAngles, int d_NoOfLevels, array& FinalPolarGrid)
{

	int debugFlag_ToPrint = 0;
	int N = d_NoOfElements - 1; /* it is always even as described in the paper*/

	array FrFT_Image_X_Cube, FrFT_Image_Y_Cube;
	array finalSeq_X, finalSeqConj_X;
	array finalSeq_Y, finalSeqConj_Y;
	for (int l = 1; l <= d_NoOfLevels; ++l)
	{
		double angle = l * Pi / d_NoOfAngles;
		array alpha = constant( cos(af::Pi / d_NoOfAngles * l),1,1 );
		array beta  = constant ( sin(af::Pi / d_NoOfAngles * l),1,1);

		//  X-axis Level
		/*-------------------- First dimension uniform FrFT for each Image per level-----------------------*/
		FrFT_Image_X_Cube = Compute2DColumnwise_FrFTUniform(Image.T(), alpha, d_NoOfElements, 1);
		FrFT_Image_X_Cube = FrFT_Image_X_Cube.T();       // Now it needs operation to the other dimension
		//af_print(FrFT_Image_X_Cube);
		
		Compute2DComplementaryLines_FrFTVariableScales((FrFT_Image_X_Cube), beta, finalSeq_X, finalSeqConj_X, d_NoOfElements, 1);
		finalSeqConj_X = flip(finalSeqConj_X, 1);             // Special operation

		/*af_print(finalSeq_X);
		af_print(finalSeqConj_X);
		af_print(FinalPolarGrid(l, span));
		af_print(FinalPolarGrid(d_NoOfAngles - l-1, span));*/
		FinalPolarGrid(l, span) = finalSeq_X;
		FinalPolarGrid(d_NoOfAngles - l, span) = finalSeqConj_X;

		if (l * 180 / d_NoOfAngles == 45)
			continue;
		
		// Y-axis level
		/*-------------------- First dimension uniform FrFT for each Image per level-----------------------*/
		FrFT_Image_Y_Cube = Compute2DColumnwise_FrFTUniform((Image), alpha, d_NoOfElements, 1);
		FrFT_Image_Y_Cube = FrFT_Image_Y_Cube.T();

		Compute2DComplementaryLines_FrFTVariableScales(FrFT_Image_Y_Cube, beta, finalSeq_Y, finalSeqConj_Y, d_NoOfElements, 1);

		FinalPolarGrid(d_NoOfAngles / 2 - l, span) = finalSeq_Y;
		FinalPolarGrid(d_NoOfAngles / 2 + l, span) = finalSeqConj_Y;
	}

	//   Computing seperately for two special indexes
	float zeroIndex = 0;
	float nintyIndex = d_NoOfAngles / 2;
	/*-------------------Precomputing the Essential Sequence Cubes :: All complex values here --------------------*/
	array lineSpacing = array(seq(-N / 2, N / 2));
	array lineSpacing_tiled_Y = tile(lineSpacing, 1, d_NoOfElements);
	array lineSpacing_Square = lineSpacing_tiled_Y * lineSpacing_tiled_Y.T();
	array ZeroNinty_Factor = cexp(complex(0, -2 * af::Pi * lineSpacing_Square * 1 / d_NoOfElements));
	if (debugFlag_ToPrint)
	{
		//af_print(lineSpacing);
		//af_print(lineSpacing_tiled_Y);
		af_print(lineSpacing_Square);
		af_print(ZeroNinty_Factor);
	}

	array ZeroLineFrFT_Image_X_Cube = FrFT_Image_Y_Cube.slice(zeroIndex).col(N / 2);
	array NintyLineFrFT_Image_Y_Cube = FrFT_Image_X_Cube.slice(zeroIndex).col(N / 2);
	if (debugFlag_ToPrint)
	{
		af_print(ZeroLineFrFT_Image_X_Cube);
		af_print(NintyLineFrFT_Image_Y_Cube);
	}

	array DFTZeroLine = sum(tile(ZeroLineFrFT_Image_X_Cube, 1, d_NoOfElements) *ZeroNinty_Factor);
	array DFTNinetyLine = sum(tile((NintyLineFrFT_Image_Y_Cube), 1, d_NoOfElements) *ZeroNinty_Factor);

	FinalPolarGrid(zeroIndex, span) = DFTZeroLine;
	FinalPolarGrid(nintyIndex, span) = DFTNinetyLine;

}

int main(int argc, char* argv[])
{
	try {
		info();
	
		//int myN = 50;
		//array testSize = randn(2 * myN, myN, c64);
		//array tiledTest = tile(testSize, 1, 1, myN / 4); 
		//array FFT_Test = fft(tiledTest);          // MAx capacity is 400
		//

		//array SimpleSeq = array(seq(1, 15));
		////af_print(SimpleSeq);
		//array FFT1DSimple = fft(SimpleSeq);
		////af_print(FFT1DSimple);
		//array IFFTSimple = ifft(FFT1DSimple);
		////af_print(IFFTSimple / pow(15, 2));

		int d_NoOfElements =  11; // 245;// 11; //9         // 2501  is the max tested
		int d_NoOfAngles = d_NoOfElements + 1;

		int d_NoOfLevels;
		if (((d_NoOfAngles - 2) % 4) != 0)
		{
			d_NoOfLevels = (d_NoOfAngles - 2) / 4 +1;
		}
		else
		{
			d_NoOfLevels = (d_NoOfAngles - 2) / 4 ;
		}
		

		array Image(d_NoOfElements, d_NoOfElements, elements_Image);

		//array Image = randu(d_NoOfElements, d_NoOfElements, f32);

		Image = Image.as(f32);
		Image = Image.T();
		array FinalPolarGrid;
		
		FullyVectorized2DPolarTransform( Image,  d_NoOfElements, d_NoOfAngles, d_NoOfLevels, FinalPolarGrid);
		printf("\nCompleted Fully Vectorized version");
		af_print(FinalPolarGrid);

		array NonVectorizedPolarGrid = constant(0, d_NoOfAngles, d_NoOfElements, c32);       // Pre allocated here
		VectorizedPerLevel2DPolarTransform(Image, d_NoOfElements, d_NoOfAngles, d_NoOfLevels, NonVectorizedPolarGrid);
		printf("\nCompleted Fully Vectorized-per-level version\n");
		af_print(NonVectorizedPolarGrid);

		array MaxAbsError = (sum(abs(FinalPolarGrid) - abs(NonVectorizedPolarGrid)));
		printf("\nMaximum Absolute error between fully vectorized and vectorized per level");
		af_print(sum(MaxAbsError));

		printf("Successfully completed the computations of 2D DFT on a full Polar Grid !\n");
		

	}
	catch (af::exception& e) {
		fprintf(stderr, "%s\n", e.what());
		throw;
	}

#ifdef WIN32 // pause in Windows
	if (!(argc == 2 && argv[1][0] == '-')) {
		printf("hit [enter]...");
		getchar();
	}
#endif
	return 0;
}
