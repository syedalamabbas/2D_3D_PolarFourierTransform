/*
* Example of how to use the mxGPUArray API in a MEX file.  This example shows
* how to write a MEX function that takes a gpuArray as input and returns a
* gpuArray output for 2D Radon solution, e.g. B=mexFunction(A).
*
* by Syed Alam Abbas, 5/25/2015
*/
#include <arrayfire.h>
#include <af/util.h>
#include "cuda_runtime.h"
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "math.h"

using namespace af;

//static const cdouble i_cdouble = { 0, 1 };
//static const  array i = constant(i_cdouble, 1, 1, c64);/* imaginary unit */
static  int isFirstRun_Uniform = 1;        // Flag to check if this is a first run on function  Compute2DColumnwise_FrFTUniform
static  int isFirstRun_Variable = 1;        // Flag to check if this is a first run on function  Compute2DComplementaryLines_FrFTVariableScales 
static  array PreMultiplicationFactor, PostMultiplicationFactor, Seq_En, BetaFactor;

enum SUPPORTED_PLATFORMS
{
	CUDA, OPENCL, CPU
};

/* Taking exp of complex numbers*/
array cexp(const array &in)
{
	if (!in.iscomplex()) return exp(in);
	return exp(real(in))*complex(cos(imag(in)), sin(imag(in)));
}

/* Multiplication of 2 complex numbers require 4 real multiplications*/
void SplitMultiplyComplex(array& A_Complex, array&  B_Complex, array& realRealPart, array& realImagPart, array& imagRealPart, array& imagImagPart)
{
	// Consider multiplication of complex numbers A_Complex = (a+ib); B_Complex = (c+id)	
	//array A_Complex, B_Complex;
	//array realRealPart;             // ac
	//array realImagPart;				// ad
	//array imagRealPart;				// bc
	//array imagImagPart;				// bd
	realRealPart = real(A_Complex)*real(B_Complex);
	realImagPart = real(A_Complex)*imag(B_Complex);
	imagRealPart = imag(A_Complex)*real(B_Complex);
	imagImagPart = imag(A_Complex)*imag(B_Complex);
}

/* Uniform FrFT for each column in Image*/
array Compute2DColumnwise_FrFTUniform(array & Image2D, array& ColumnScales_1D, int& d_NoOfElements_, int& d_NoOfScales)
{
	/*-----------------------------------Preparing Padded & Tiled Imag2D --------------------------------------------*/
	array Zeros = constant(0, d_NoOfElements_, d_NoOfElements_, f64);            // Generates on the device
	array Zero_Padded_Image2D = join(0, Image2D, Zeros);
	array Image2D_Tiled = tile(Zero_Padded_Image2D, 1, 1, d_NoOfScales);
	int N = d_NoOfElements_ - 1;
	if (isFirstRun_Uniform == 1)
	{
		/*-------------------------------------------Creating Index Cubes and Sequences----------------------------------------------------*/
		array leftSideIndexes = array(seq(0, N, 1)).as(f64);
		array rightSideIndexesOnes = -1 * array(seq(1, d_NoOfElements_, 1)).as(f64);
		array rightSideIndexesZeros = constant(0, d_NoOfElements_, 1, f64);
		array rightSideIndexesN_2 = constant(N / 2, d_NoOfElements_, 1, f64);

		array indexedElementsEn = join(0, leftSideIndexes, flip(rightSideIndexesOnes, 0));
		array indexedElementsPre = join(0, leftSideIndexes, rightSideIndexesN_2);     /* This is for Keeping pre and post multiplication factor upper half only*/
		array indexedElementsPost = join(0, leftSideIndexes, rightSideIndexesZeros);

		array indexedElements_Tiled_En = tile(pow(indexedElementsEn, 2), 1, d_NoOfElements_, d_NoOfScales);
		array indexedElements_Tiled_PreMulti = tile(indexedElementsPre - N / 2, 1, d_NoOfElements_, d_NoOfScales);
		array indexedElements_Tiled_PostMulti = tile(indexedElementsPost, 1, d_NoOfElements_, d_NoOfScales);

		/*--------------------------Creating FrFT scale cubes------------------------------------*/
		array ColumnScales_1D_Mods = moddims(ColumnScales_1D, 1, 1, d_NoOfScales);
		array ColumnScales_1D_Tiled_depth = tile(ColumnScales_1D_Mods, 2 * d_NoOfElements_, d_NoOfElements_, 1);
		array ColumnScales_1D_FullTiled = moddims(ColumnScales_1D_Tiled_depth, 2 * d_NoOfElements_, d_NoOfElements_, d_NoOfScales);

		/*-------------------Precomputing the Essential Sequence Cubes :: All complex values here --------------------*/
		//array imaginaryUnit_Tiled = tile(i, 2 * d_NoOfElements_, d_NoOfElements_, d_NoOfScales);
		Seq_En = cexp( complex(0, - af::Pi * indexedElements_Tiled_En * ColumnScales_1D_FullTiled / d_NoOfElements_));   /* E(n) as defined in the paper*/
		array Ones = constant(1, d_NoOfElements_, d_NoOfElements_, f64);
		array subtractValues = tile(join(0, Zeros, Ones), 1, 1, d_NoOfScales);			/* This is for  Keeping pre and post multiplication factor upper half only*/
		PreMultiplicationFactor = cexp( complex (0, af::Pi * indexedElements_Tiled_PreMulti * ColumnScales_1D_FullTiled * N / d_NoOfElements_)) - subtractValues;
		PostMultiplicationFactor = cexp( complex ( 0,  af::Pi *  indexedElements_Tiled_PostMulti * ColumnScales_1D_FullTiled * N / d_NoOfElements_)) - subtractValues;
		isFirstRun_Uniform = 0;
		//af::deviceGC();
	}
	/*--------------------Preprocessing Cubes-----------------------*/
	array Image2D_Tiled_PreMulti = Image2D_Tiled * PreMultiplicationFactor;
	array Image2D_Tiled_PreMulti_SeqEn = Image2D_Tiled_PreMulti * Seq_En;


	/*-------------------Computing Convolution--------------------*/
	array firstFFT_X = fft(Image2D_Tiled_PreMulti_SeqEn);
	array secondFFT_X = fft(conjg(Seq_En));
	array interim_FrFT_X = ifft(firstFFT_X * secondFFT_X);


	/*-------------------Postprocessing-----------------------------*/
	array  FrFT_Image_X = interim_FrFT_X * Seq_En * PostMultiplicationFactor;

	/*--------------------Grab only the top half elements drop overlapping------------------*/
	array FrFT_Image_X_Cube = FrFT_Image_X.rows(0, N);
	return FrFT_Image_X_Cube;
}

/* Variable FrFT for each column in Image*/
void Compute2DComplementaryLines_FrFTVariableScales(array & OneD_FrFT, array& ColumnScales_1D, array& final2DFrFTImage, array& final2DFrFTConjImage, int d_NoOfElements, int d_NoOfScales)
{
	int N = d_NoOfElements - 1;
	if (isFirstRun_Variable == 1)
	{
		array lineSpacing = array(seq(-N / 2, N / 2)).as(f64);
		array lineSpacing_tiled_Y = tile(lineSpacing, 1, d_NoOfElements);
		array lineSpacing_Square = lineSpacing_tiled_Y * lineSpacing_tiled_Y.T();
		array lineSpacing_Square_TiledLevel = tile(lineSpacing_Square, 1, 1, d_NoOfScales);
	
		//af_print(beta_Levels);
		array beta_Mods = moddims(ColumnScales_1D, 1, 1, d_NoOfScales);
		array beta_Tiled_depth = tile(beta_Mods, d_NoOfElements, d_NoOfElements, 1);
		//af_print(beta_Tiled_depth);
		array beta_Tiled = moddims(beta_Tiled_depth, d_NoOfElements, d_NoOfElements, d_NoOfScales);
	
		/*-------------------Precomputing the Essential Sequence Cubes :: All complex values here --------------------*/
		BetaFactor = cexp( complex(0, -2 * af::Pi * lineSpacing_Square_TiledLevel * beta_Tiled / d_NoOfElements));
	
		isFirstRun_Variable = 0;
		af::deviceGC();
	}
	// Consider multiplication of complex numbers A = (a+ib); B = (c+id)	
	array realRealPart;             // ac
	array realImagPart;				// ad 
	array imagRealPart;				// bc
	array imagImagPart;				// bd

	SplitMultiplyComplex(OneD_FrFT, BetaFactor, realRealPart, realImagPart, imagRealPart, imagImagPart);

	array tempSeq_X = sum(complex(realRealPart - imagImagPart, realImagPart + imagRealPart)); // sum(flip(FrFT_Image_X_Cube, 0) *BetaFactor);
	array tempSeqConj_X = sum(complex(realRealPart + imagImagPart, imagRealPart - realImagPart));//  sum(flip(FrFT_Image_X_Cube, 0) * conjg(BetaFactor));
	
	final2DFrFTImage = moddims(tempSeq_X, d_NoOfElements, d_NoOfScales).T();
	final2DFrFTConjImage = moddims(tempSeqConj_X, d_NoOfElements, d_NoOfScales).T();
}

/*
* High level Host code
* Computes the FrFT centered using the definition given in the paper,
"An exact and fast computation of Discrete Fourier Transform for polar grid and spherical grid"
by Syed Alam Abbas, 5/25/2015
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
	try {
		/* Initialize the MathWorks GPU API. */ 
		mxInitGPU();

		mexPrintf("Executing custom mex for computing 2D DFT on a Polar Grid using ArrayFire GPU accelerated library latest!");

		// Validate the input
		if (nrhs < 5 || nlhs < 2) {
			mexErrMsgTxt("Expected 5 inputs and 2 output.");
		}

		/*Input Variables*/
		const double* d_Image;
		int d_NoOfAngles;
		int d_NoOfLevels;
		int d_NoOfElements;
		mxGPUArray const * MxInputImage;
		mxArray* MxInputImageCPU;
		
		int PLATFORM = (size_t)mxGetScalar(prhs[4]);        // Given as an input PLATFORM

		switch (PLATFORM)                   // The settings change for input
		{
		case CUDA:
		case OPENCL:
			MxInputImage = mxGPUCreateFromMxArray(prhs[0]);   // GPU
			/* extract a pointer to the input data which is a real image on the device.*/
			d_Image = (double const *)(mxGPUGetDataReadOnly(MxInputImage));  // GPU 
			break;
		case CPU:
			MxInputImageCPU = mxDuplicateArray(prhs[0]);
			d_Image = mxGetPr(MxInputImageCPU);
			break;
		default:
			break;
		}


		/* Collect the input data from MATLAB MexArray RHS */
		d_NoOfAngles = (size_t)mxGetScalar(prhs[1]);            /*Check it, this should always be even*/
		d_NoOfLevels = (size_t)mxGetScalar(prhs[2]);
		d_NoOfElements = (size_t)mxGetScalar(prhs[3]);			/*Check it, this should always be odd*/
		int N = d_NoOfElements - 1;								/* it is always even as described in the paper*/

		/*********************Creating Array Fire object************************************/
		array Image(d_NoOfElements, d_NoOfElements, d_Image);


		/*--------------------------Creating Alpha levels------------------------------------*/
		array alpha_Levels = cos(af::Pi / ((double) d_NoOfAngles) * array(seq(1, d_NoOfLevels)).as(f64));

		/*--------------------------Creating Beta levels---------------------------------------*/
		array beta_Levels = sin(Pi / ((double)d_NoOfAngles) * array(seq(1, d_NoOfLevels)).as(f64));

		/*-------------------Precomputing the Essential Sequence Cubes :: All complex values here --------------------*/
		array lineSpacing = array(seq(-N / 2, N / 2)).as(f64);
		array lineSpacing_tiled_Y = tile(lineSpacing, 1, d_NoOfElements);
		array lineSpacing_Square = lineSpacing_tiled_Y * lineSpacing_tiled_Y.T();
		array ZeroNinty_Factor = cexp( complex( 0, -2 * af::Pi * lineSpacing_Square * 1 / d_NoOfElements));

		/*-------------------- First dimension uniform FrFT for each Image per level-----------------------*/
		array FrFT_Image_X_Cube = Compute2DColumnwise_FrFTUniform(Image.T(), alpha_Levels, d_NoOfElements, d_NoOfLevels);
		switch (PLATFORM)
		{
		case CUDA:
			af::deviceGC();
		default:
			break; 
		}

		array FrFT_Image_Y_Cube = Compute2DColumnwise_FrFTUniform((Image), alpha_Levels, d_NoOfElements, d_NoOfLevels);

		FrFT_Image_X_Cube = FrFT_Image_X_Cube.T();       // Now it needs operation to the other dimension
		FrFT_Image_Y_Cube = FrFT_Image_Y_Cube.T();

		
		/*--------------------Finally all computations for  the  Polar Grid-----------*/
		//   Computing for all the grid expect two special indexes
		array levelSeq = array(seq(0, d_NoOfLevels - 1)).as(f64);
		array finalIndexSeq1_X = 1 + levelSeq;
		array finalIndexSeq2_X = d_NoOfAngles - finalIndexSeq1_X;
		array finalIndexSeq3_Y = d_NoOfAngles / 2 - finalIndexSeq1_X;
		array finalIndexSeq4_Y = d_NoOfAngles / 2 + finalIndexSeq1_X;


		array finalSeq_X, finalSeqConj_X;
		Compute2DComplementaryLines_FrFTVariableScales((FrFT_Image_X_Cube), beta_Levels, finalSeq_X, finalSeqConj_X, d_NoOfElements, d_NoOfLevels);
		finalSeqConj_X = flip(finalSeqConj_X, 1);             // Special operation
		array finalSeq_Y, finalSeqConj_Y;
		Compute2DComplementaryLines_FrFTVariableScales(FrFT_Image_Y_Cube, beta_Levels, finalSeq_Y, finalSeqConj_Y, d_NoOfElements, d_NoOfLevels);

		// Removing just 2 redundant computations for 45 degree case
		if (0 == remainder(d_NoOfAngles, 4))
		{
			finalIndexSeq3_Y = finalIndexSeq3_Y.rows(0, d_NoOfLevels - 2);          // Removing just the last rows from 4 structures
			finalSeq_Y = finalSeq_Y.rows(0, d_NoOfLevels - 2);
			finalIndexSeq4_Y = finalIndexSeq4_Y.rows(0, d_NoOfLevels - 2);
			finalSeqConj_Y = finalSeqConj_Y.rows(0, d_NoOfLevels - 2);
		}

		//   Computing seperately for two special indexes
		double zeroIndex = 0;
		double nintyIndex = d_NoOfAngles / 2;
		double values[] = { zeroIndex, nintyIndex };
		array SpecialTwoIndexes(2, 1, values);

		array ZeroLineFrFT_Image_X_Cube = FrFT_Image_Y_Cube.slice(zeroIndex).col(N / 2);
		array NintyLineFrFT_Image_Y_Cube = FrFT_Image_X_Cube.slice(zeroIndex).col(N / 2);

		array DFTZeroLine = sum(tile(ZeroLineFrFT_Image_X_Cube, 1, d_NoOfElements) *ZeroNinty_Factor);
		array DFTNinetyLine = sum(tile((NintyLineFrFT_Image_Y_Cube), 1, d_NoOfElements) *ZeroNinty_Factor);
		array SpecialTwoLines = join(0, DFTZeroLine, DFTNinetyLine);

		array UnsortedIndexes = join(0, join(0, join(0, join(0, finalIndexSeq1_X, finalIndexSeq2_X), finalIndexSeq3_Y), finalIndexSeq4_Y), SpecialTwoIndexes);
		array tiledUnsortedIndexes = tile(UnsortedIndexes, 1, d_NoOfElements);
		array UnsortedPolarGrid = join(0, join(0, join(0, join(0, finalSeq_X, finalSeqConj_X), finalSeq_Y), finalSeqConj_Y), SpecialTwoLines);


		array FinalPolarGridReal;// = constant(0, d_NoOfElements, d_NoOfAngles, c64);
		array Output_Keys_Sorted;
		sort(Output_Keys_Sorted, FinalPolarGridReal, tiledUnsortedIndexes, real(UnsortedPolarGrid));

		array FinalPolarGridImag;// = constant(0, d_NoOfElements, d_NoOfAngles, c64);
		array Output_Keys_Sorted2;
		sort(Output_Keys_Sorted2, FinalPolarGridImag, tiledUnsortedIndexes, imag(UnsortedPolarGrid));


		mexPrintf("\nSuccessfully completed the computations of 2D DFT on a full Polar Grid  %d-by-%d!", d_NoOfAngles,d_NoOfElements);
		
		double* d_FinalPolarGridReal;       // Device pointer obtained from ArrayFire computations
		double* d_FinalPolarGridImag;		// Device pointer obtained from ArrayFire computations
		double* PolarGridReal_OUTPUT;       // MATLAB output pointer to be copied to the solution
		double* PolarGridImag_OUTPUT;		// MATLAB output pointer to be copied to the solution

		mwSize dims[] = { d_NoOfAngles, d_NoOfElements };

		switch (PLATFORM)                   // The settings change for input
		{
		case CUDA:
		case OPENCL:
			// Final processed double pointers
			d_FinalPolarGridReal = FinalPolarGridReal.device<double>();
			d_FinalPolarGridImag = FinalPolarGridImag.device<double>();

			/*output variables*/
			mxGPUArray* mxOutputRealPolarGridImage;
			mxGPUArray* mxOutputImagPolarGridImage;
			

			/* Create a GPUArray to hold the result and get its underlying pointer. */
			mxOutputRealPolarGridImage = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(MxInputImage),
				dims,
				mxGPUGetClassID(MxInputImage),
				mxGPUGetComplexity(MxInputImage),
				MX_GPU_DO_NOT_INITIALIZE);
			PolarGridReal_OUTPUT = (double *)(mxGPUGetData(mxOutputRealPolarGridImage));

			mxOutputImagPolarGridImage = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(MxInputImage),
				dims,
				mxGPUGetClassID(MxInputImage),
				mxGPUGetComplexity(MxInputImage),
				MX_GPU_DO_NOT_INITIALIZE);
			PolarGridImag_OUTPUT = (double *)(mxGPUGetData(mxOutputImagPolarGridImage));
			/* Copy processed Values from array object to MxArrayRealData */
			cudaMemcpy(PolarGridReal_OUTPUT, d_FinalPolarGridReal, d_NoOfAngles*d_NoOfElements* sizeof(double), cudaMemcpyDeviceToDevice);
			cudaMemcpy(PolarGridImag_OUTPUT, d_FinalPolarGridImag, d_NoOfAngles*d_NoOfElements* sizeof(double), cudaMemcpyDeviceToDevice);

			/* Wrap the result up as a MATLAB gpuArray for return. */
			plhs[0] = mxGPUCreateMxArrayOnGPU(mxOutputRealPolarGridImage);
			plhs[1] = mxGPUCreateMxArrayOnGPU(mxOutputImagPolarGridImage);
			/*
			* The mxGPUArray pointers are host-side structures that refer to device
			* data. These must be destroyed before leaving the MEX function.
			*/
			mxGPUDestroyGPUArray(MxInputImage);

			break;
		case CPU:
			// Final processed double pointers
			d_FinalPolarGridReal = FinalPolarGridReal.host<double>();   // Source
			d_FinalPolarGridImag = FinalPolarGridImag.host<double>();

			mxArray*  mxOutputRealPolarGridImageCPU;
			mxArray*  mxOutputImagPolarGridImageCPU;
			mxOutputRealPolarGridImageCPU = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
			mxOutputImagPolarGridImageCPU = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
			PolarGridReal_OUTPUT = mxGetPr(mxOutputRealPolarGridImageCPU);
			PolarGridImag_OUTPUT = mxGetPr(mxOutputImagPolarGridImageCPU);

			memcpy(PolarGridReal_OUTPUT, d_FinalPolarGridReal, d_NoOfAngles*d_NoOfElements* sizeof(double));
			memcpy(PolarGridImag_OUTPUT, d_FinalPolarGridImag, d_NoOfAngles*d_NoOfElements* sizeof(double));

			plhs[0] = mxOutputRealPolarGridImageCPU;
			plhs[1] = mxOutputImagPolarGridImageCPU;
			break;
		default:
			break;
		}


		mexPrintf("\nFinished processing custom CUDA mex with ArrayFire for computing 2D DFT on Polar Grid, Status = Success\n");

	}
	catch (af::exception &ex) {
		mexPrintf("%s\n", ex.what());
	}

	
}
 