/*
* This is a program for testing the 3D Spherical Transform
*/

#include <arrayfire.h>
#include <af/util.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <windows.h>
#include <tlhelp32.h>
#include <stdio.h>

using namespace af;


enum AXES
{
	X_AXIS = 1,
	Y_AXIS = 0,   // This is how arrayFire refers to Y-axis
	Z_AXIS = 2
};

static const af::dtype PRECISION_REAL = f64;                   // Change here from single to double
static const af::dtype PRECISION_COMPLEX = c64;


array cexp(const array &in)
{
	if (!in.iscomplex()) return exp(in);
	return exp(real(in))*complex(cos(imag(in)), sin(imag(in)));
}

// This is a test image constant so we can compare exact values of the lines data computed

static const double ThreeDImageValues[] = 
{
	0.8147, 0.9649, 0.7922, 0.3922, 0.6948, 0.4898, 0.1190, 0.6991, 0.8143,
	0.9058, 0.1576, 0.9595, 0.6555, 0.3171, 0.4456, 0.4984, 0.8909, 0.2435,
	0.1270, 0.9706, 0.6557, 0.1712, 0.9502, 0.6463, 0.9597, 0.9593, 0.9293,
	0.9134, 0.9572, 0.0357, 0.7060, 0.0344, 0.7094, 0.3404, 0.5472, 0.3500,
	0.6324, 0.4854, 0.8491, 0.0318, 0.4387, 0.7547, 0.5853, 0.1386, 0.1966,
	0.0975, 0.8003, 0.9340, 0.2769, 0.3816, 0.2760, 0.2238, 0.1493, 0.2511,
	0.2785, 0.1419, 0.6787, 0.0462, 0.7655, 0.6797, 0.7513, 0.2575, 0.6160,
	0.5469, 0.4218, 0.7577, 0.0971, 0.7952, 0.6551, 0.2551, 0.8407, 0.4733,
	0.9575, 0.9157, 0.7431, 0.8235, 0.1869, 0.1626, 0.5060, 0.2543, 0.3517, 

	//I(:, : , 2) =

	0.8308, 0.0759, 0.3371, 0.6892, 0.9961, 0.0844, 0.1361, 0.4018, 0.9448,
	0.5853, 0.0540, 0.1622, 0.7482, 0.0782, 0.3998, 0.8693, 0.0760, 0.4909,
	0.5497, 0.5308, 0.7943, 0.4505, 0.4427, 0.2599, 0.5797, 0.2399, 0.4893,
	0.9172, 0.7792, 0.3112, 0.0838, 0.1067, 0.8001, 0.5499, 0.1233, 0.3377,
	0.2858, 0.9340, 0.5285, 0.2290, 0.9619, 0.4314, 0.1450, 0.1839, 0.9001,
	0.7572, 0.1299, 0.1656, 0.9133, 0.0046, 0.9106, 0.8530, 0.2400, 0.3692,
	0.7537, 0.5688, 0.6020, 0.1524, 0.7749, 0.1818, 0.6221, 0.4173, 0.1112,
	0.3804, 0.4694, 0.2630, 0.8258, 0.8173, 0.2638, 0.3510, 0.0497, 0.7803,
	0.5678, 0.0119, 0.6541, 0.5383, 0.8687, 0.1455, 0.5132, 0.9027, 0.3897,

	//I(:, : , 3) =

	0.2417, 0.3532, 0.5470, 0.0811, 0.8176, 0.5502, 0.2259, 0.9797, 0.2217,
	0.4039, 0.8212, 0.2963, 0.9294, 0.7948, 0.6225, 0.1707, 0.4389, 0.1174,
	0.0965, 0.0154, 0.7447, 0.7757, 0.6443, 0.5870, 0.2277, 0.1111, 0.2967,
	0.1320, 0.0430, 0.1890, 0.4868, 0.3786, 0.2077, 0.4357, 0.2581, 0.3188,
	0.9421, 0.1690, 0.6868, 0.4359, 0.8116, 0.3012, 0.3111, 0.4087, 0.4242,
	0.9561, 0.6491, 0.1835, 0.4468, 0.5328, 0.4709, 0.9234, 0.5949, 0.5079,
	0.5752, 0.7317, 0.3685, 0.3063, 0.3507, 0.2305, 0.4302, 0.2622, 0.0855,
	0.0598, 0.6477, 0.6256, 0.5085, 0.9390, 0.8443, 0.1848, 0.6028, 0.2625,
	0.2348, 0.4509, 0.7802, 0.5108, 0.8759, 0.1948, 0.9049, 0.7112, 0.8010,

	//I(:, : , 4) =

	0.0292, 0.5211, 0.8852, 0.1068, 0.1978, 0.8055, 0.7127, 0.8181, 0.4538,
	0.9289, 0.2316, 0.9133, 0.6538, 0.0305, 0.5767, 0.5005, 0.8175, 0.4324,
	0.7303, 0.4889, 0.7962, 0.4942, 0.7441, 0.1829, 0.4711, 0.7224, 0.8253,
	0.4886, 0.6241, 0.0987, 0.7791, 0.5000, 0.2399, 0.0596, 0.1499, 0.0835,
	0.5785, 0.6791, 0.2619, 0.7150, 0.4799, 0.8865, 0.6820, 0.6596, 0.1332,
	0.2373, 0.3955, 0.3354, 0.9037, 0.9047, 0.0287, 0.0424, 0.5186, 0.1734,
	0.4588, 0.3674, 0.6797, 0.8909, 0.6099, 0.4899, 0.0714, 0.9730, 0.3909,
	0.9631, 0.9880, 0.1366, 0.3342, 0.6177, 0.1679, 0.5216, 0.6490, 0.8314,
	0.5468, 0.0377, 0.7212, 0.6987, 0.8594, 0.9787, 0.0967, 0.8003, 0.8034,

	//I(:, : , 5) =

	0.0605, 0.9841, 0.0527, 0.7011, 0.0326, 0.8555, 0.3846, 0.3439, 0.4253,
	0.3993, 0.1672, 0.7379, 0.6663, 0.5612, 0.6448, 0.5830, 0.5841, 0.3127,
	0.5269, 0.1062, 0.2691, 0.5391, 0.8819, 0.3763, 0.2518, 0.1078, 0.1615,
	0.4168, 0.3724, 0.4228, 0.6981, 0.6692, 0.1909, 0.2904, 0.9063, 0.1788,
	0.6569, 0.1981, 0.5479, 0.6665, 0.1904, 0.4283, 0.6171, 0.8797, 0.4229,
	0.6280, 0.4897, 0.9427, 0.1781, 0.3689, 0.4820, 0.2653, 0.8178, 0.0942,
	0.2920, 0.3395, 0.4177, 0.1280, 0.4607, 0.1206, 0.8244, 0.2607, 0.5985,
	0.4317, 0.9516, 0.9831, 0.9991, 0.9816, 0.5895, 0.9827, 0.5944, 0.4709,
	0.0155, 0.9203, 0.3015, 0.1711, 0.1564, 0.2262, 0.7302, 0.0225, 0.6959,

	//I(:, : , 6) =

	0.6999, 0.7184, 0.2665, 0.6377, 0.2240, 0.9160, 0.0358, 0.2428, 0.5466,
	0.6385, 0.9686, 0.1537, 0.9577, 0.6678, 0.0012, 0.1759, 0.9174, 0.4257,
	0.0336, 0.5313, 0.2810, 0.2407, 0.8444, 0.4624, 0.7218, 0.2691, 0.6444,
	0.0688, 0.3251, 0.4401, 0.6761, 0.3445, 0.4243, 0.4735, 0.7655, 0.6476,
	0.3196, 0.1056, 0.5271, 0.2891, 0.7805, 0.4609, 0.1527, 0.1887, 0.6790,
	0.5309, 0.6110, 0.4574, 0.6718, 0.6753, 0.7702, 0.3411, 0.2875, 0.6358,
	0.6544, 0.7788, 0.8754, 0.6951, 0.0067, 0.3225, 0.6074, 0.0911, 0.9452,
	0.4076, 0.4235, 0.5181, 0.0680, 0.6022, 0.7847, 0.1917, 0.5762, 0.2089,
	0.8200, 0.0908, 0.9436, 0.2548, 0.3868, 0.4714, 0.7384, 0.6834, 0.7093,

	//I(:, : , 7) =

	0.2362, 0.4162, 0.3181, 0.7210, 0.3658, 0.0938, 0.3477, 0.3592, 0.2703,
	0.1194, 0.8419, 0.1192, 0.5225, 0.7635, 0.5254, 0.1500, 0.7363, 0.1971,
	0.6073, 0.8329, 0.9398, 0.9937, 0.6279, 0.5303, 0.5861, 0.3947, 0.8217,
	0.4501, 0.2564, 0.6456, 0.2187, 0.7720, 0.8611, 0.2621, 0.6834, 0.4299,
	0.4587, 0.6135, 0.4795, 0.1058, 0.9329, 0.4849, 0.0445, 0.7040, 0.8878,
	0.6619, 0.5822, 0.6393, 0.1097, 0.9727, 0.3935, 0.7549, 0.4423, 0.3912,
	0.7703, 0.5407, 0.5447, 0.0636, 0.1920, 0.6714, 0.2428, 0.0196, 0.7691,
	0.3502, 0.8699, 0.6473, 0.4046, 0.1389, 0.7413, 0.4424, 0.3309, 0.3968,
	0.6620, 0.2648, 0.5439, 0.4484, 0.6963, 0.5201, 0.6878, 0.4243, 0.8085,

	//I(:, : , 8) =

	0.7551, 0.7689, 0.4070, 0.6787, 0.6967, 0.5277, 0.5860, 0.7690, 0.2094,
	0.3774, 0.1673, 0.7487, 0.4952, 0.5828, 0.4795, 0.2467, 0.5814, 0.5523,
	0.2160, 0.8620, 0.8256, 0.1897, 0.8154, 0.8013, 0.6664, 0.9283, 0.6299,
	0.7904, 0.9899, 0.7900, 0.4950, 0.8790, 0.2278, 0.0835, 0.5801, 0.0320,
	0.9493, 0.5144, 0.3185, 0.1476, 0.9889, 0.4981, 0.6260, 0.0170, 0.6147,
	0.3276, 0.8843, 0.5341, 0.0550, 0.0005, 0.9009, 0.6609, 0.1209, 0.3624,
	0.6713, 0.5880, 0.0900, 0.8507, 0.8654, 0.5747, 0.7298, 0.8627, 0.0495,
	0.4386, 0.1548, 0.1117, 0.5606, 0.6126, 0.8452, 0.8908, 0.4843, 0.4896,
	0.8335, 0.1999, 0.1363, 0.9296, 0.9900, 0.7386, 0.9823, 0.8449, 0.1925,

	//I(:, : , 9) =

	0.1231, 0.4991, 0.5650, 0.6210, 0.9844, 0.3013, 0.3479, 0.5400, 0.1781,
	0.2055, 0.5358, 0.6403, 0.5737, 0.8589, 0.2955, 0.4460, 0.7069, 0.3596,
	0.1465, 0.4452, 0.4170, 0.0521, 0.7856, 0.3329, 0.0542, 0.9995, 0.0567,
	0.1891, 0.1239, 0.2060, 0.9312, 0.5134, 0.4671, 0.1771, 0.2878, 0.5219,
	0.0427, 0.4904, 0.9479, 0.7287, 0.1776, 0.6482, 0.6628, 0.4145, 0.3358,
	0.6352, 0.8530, 0.0821, 0.7378, 0.3986, 0.0252, 0.3308, 0.4648, 0.1757,
	0.2819, 0.8739, 0.1057, 0.0634, 0.1339, 0.8422, 0.8985, 0.7640, 0.2089,
	0.5386, 0.2703, 0.1420, 0.8604, 0.0309, 0.5590, 0.1182, 0.8182, 0.9052,
	0.6952, 0.2085, 0.1665, 0.9344, 0.9391, 0.8541, 0.9884, 0.1002, 0.6754
};




void SplitMultiplyComplex(array& A, array&  B, array& realRealPart, array& realImagPart, array& imagRealPart, array& imagImagPart)
{
	// Consider multiplication of complex numbers A = (a+ib); B = (c+id)	
	//array A, B;
	//array realRealPart;             // ac
	//array realImagPart;				// ad
	//array imagRealPart;				// bc
	//array imagImagPart;				// bd
	realRealPart = real(A)*real(B);
	//af_print(realRealPart);
	realImagPart = real(A)*imag(B);
	imagRealPart = imag(A)*real(B);
	imagImagPart = imag(A)*imag(B);
}


static array indexedElements_Tiled_En;
static array indexedElements_Tiled_PreMulti;
static array indexedElements_Tiled_PostMulti;

void GlobalArraysComputeInitialize(int d_NoOfElements)
{
	/*-------------------------------------------Creating Index Cubes and Sequences----------------------------------------------------*/
	int N = d_NoOfElements - 1;
	array leftSideIndexes = array(seq(0, N, 1)).as(PRECISION_REAL);
	// af_print(leftSideIndexes);
	array rightSideIndexesOnes = -1 * array(seq(1, d_NoOfElements, 1)).as(PRECISION_REAL);
	//af_print(flip( rightSideIndexesOnes,0));
	array rightSideIndexesZeros = constant(0, d_NoOfElements, 1, PRECISION_REAL);
	//af_print(rightSideIndexesZeros);
	array rightSideIndexesN_2 = constant(N / 2, d_NoOfElements, 1, PRECISION_REAL);
	//af_print(rightSideIndexesN_2);

	array indexedElementsEn = join(0, leftSideIndexes, flip(rightSideIndexesOnes, 0));
	array indexedElementsPre = join(0, leftSideIndexes, rightSideIndexesN_2);     /* This is for Keeping pre and post multiplication factor upper half only*/
	array indexedElementsPost = join(0, leftSideIndexes, rightSideIndexesZeros);
	//af_print(indexedElementsEn);
	//af_print(indexedElementsPre_Post);

	 indexedElements_Tiled_En = tile(pow(indexedElementsEn, 2), 1, d_NoOfElements, d_NoOfElements);
	//af_print(indexedElements_Tiled_En.slice(0));
	 indexedElements_Tiled_PreMulti = tile(indexedElementsPre - N / 2, 1, d_NoOfElements, d_NoOfElements);
	//af_print(indexedElements_Tiled_PreMulti.slice(0));
	 indexedElements_Tiled_PostMulti = tile(indexedElementsPost, 1, d_NoOfElements, d_NoOfElements);
	//af_print(indexedElements_Tiled_PostMulti.slice(0));
}


array ComputeNonVectorColumnwise_FrFTUniform(array& Image_3D, double ColumnScale, int d_NoOfElements)
{
	// Image_3D is unpadded just a 3D volume (d_NoOfElements)^3 whose column wise FrFT scaling needs to be done
	// ColumnScale is the given scale that needs to be operated on the image uniformly for entire volume (3D) 
	// d_NoOfElements is the number of elements in the volume Image_3D which is a perfect cube, 


	/**************************************************Preparing the 3D Volume Stack***************************************************/
	array zeros = constant(0, d_NoOfElements, d_NoOfElements, d_NoOfElements, PRECISION_REAL);
	array paddedImage_3D = join(0, Image_3D, zeros);									  // Column Wise padding
	//af_print(paddedImage_3D.slice(0));
	int N = d_NoOfElements - 1;

	/*-------------------------------------------Creating Index Cubes and Sequences----------------------------------------------------*/
	
	//array leftSideIndexes = array(seq(0, N, 1)).as(PRECISION_REAL);
	//// af_print(leftSideIndexes);
	//array rightSideIndexesOnes = -1 * array(seq(1, d_NoOfElements, 1)).as(PRECISION_REAL);
	////af_print(flip( rightSideIndexesOnes,0));
	//array rightSideIndexesZeros = constant(0, d_NoOfElements, 1, PRECISION_REAL);
	////af_print(rightSideIndexesZeros);
	//array rightSideIndexesN_2 = constant(N / 2, d_NoOfElements, 1, PRECISION_REAL);
	////af_print(rightSideIndexesN_2);

	//array indexedElementsEn = join(0, leftSideIndexes, flip(rightSideIndexesOnes, 0));
	//array indexedElementsPre = join(0, leftSideIndexes, rightSideIndexesN_2);     /* This is for Keeping pre and post multiplication factor upper half only*/
	//array indexedElementsPost = join(0, leftSideIndexes, rightSideIndexesZeros);
	////af_print(indexedElementsEn);
	////af_print(indexedElementsPre_Post);

	//array indexedElements_Tiled_En = tile(pow(indexedElementsEn, 2), 1, d_NoOfElements, d_NoOfElements);
	////af_print(indexedElements_Tiled_En.slice(0));
	//array indexedElements_Tiled_PreMulti = tile(indexedElementsPre - N / 2, 1, d_NoOfElements, d_NoOfElements);
	////af_print(indexedElements_Tiled_PreMulti.slice(0));
	//array indexedElements_Tiled_PostMulti = tile(indexedElementsPost, 1, d_NoOfElements, d_NoOfElements);
	////af_print(indexedElements_Tiled_PostMulti.slice(0));


	/*-------------------Precomputing the Essential Sequence Cubes :: All complex values here --------------------*/
	//array ImaginaryUnitTiled = tile(i, 2 * d_NoOfElements, d_NoOfElements, d_NoOfElements);
	//af_print(ImaginaryUnitTiled.slice(0));
	array Seq_En = cexp(complex(0, -Pi * indexedElements_Tiled_En * ColumnScale / d_NoOfElements));   /* E(n) as defined in the paper*/
	//af_print((Pi * indexedElements_Tiled_En * ColumnScale / d_NoOfElements).slice(0));
	//af_print(Seq_En.slice(0));
	//array ones = constant(1, d_NoOfElements, 1, PRECISION_REAL);
	//array subtractValues = tile(join(0, rightSideIndexesZeros, ones), 1, d_NoOfElements, d_NoOfElements);			/* This is for  Keeping pre and post multiplication factor upper half only*/
	array PreMultiplicationFactor = cexp(complex(0, af::Pi * indexedElements_Tiled_PreMulti * ColumnScale * N / d_NoOfElements));// -subtractValues;
	//af_print(PreMultiplicationFactor.slice(0));
	array PostMultiplicationFactor = cexp(complex(0, af::Pi *  indexedElements_Tiled_PostMulti * ColumnScale * N / d_NoOfElements));// - subtractValues;
	//af_print(PostMultiplicationFactor.slice(0));

	/*--------------------Preprocessing Cubes-----------------------*/
	array Image_PreMulti = paddedImage_3D * PreMultiplicationFactor;
	//af_print(Image_PreMulti.slice(0));
	array Image_PreMulti_SeqEn = Image_PreMulti * Seq_En;
	//af_print(Image_PreMulti_SeqEn.slice(0));

	/*-------------------Computing Convolution--------------------*/
	array firstFFT_Image = fft(Image_PreMulti_SeqEn);
	//af_print(firstFFT_Image.slice(0));
	array secondFFT_Seq_En = fft(conjg(Seq_En));
	//af_print(secondFFT_Seq_En.slice(0));
	array interim_FrFT_Image = ifft(firstFFT_Image * secondFFT_Seq_En);
	//af_print(interim_FrFT_Image.slice(0));

	/*-------------------Postprocessing-----------------------------*/
	array  Final_FrFT_Image = interim_FrFT_Image * Seq_En * PostMultiplicationFactor;
	//af_print(Final_FrFT_Image.slice(0));

	/*--------------------Grab only the top half elements drop overlapping------------------*/
	Final_FrFT_Image = Final_FrFT_Image.rows(0, N);
	//af_print(Final_FrFT_Image.slice(0));

	return Final_FrFT_Image;
}

void ComputeNonVector_NonUniformFrFT_SingleColumn(array&FrFT1D_Image_3D, double ColumnScale, int d_NoOfElements, array& Slice2D, array& Slice2D_Conj)
{
	int N = d_NoOfElements - 1;
	array lineSpacing = array(seq(-N / 2, N / 2)).as(PRECISION_REAL);
	array lineSpacing_tiled_Y = tile(lineSpacing, 1, d_NoOfElements);
	//af_print(lineSpacing_tiled_Y);
	array lineSpacing_Square = lineSpacing_tiled_Y * lineSpacing_tiled_Y.T();
	//af_print(lineSpacing_Square);
	array lineSpacing_Square_TiledLevel = tile(lineSpacing_Square, 1, 1, d_NoOfElements);
	//af_print(lineSpacing_Square_TiledLevel);
	array Column_ScaleFactor = cexp(complex(0, -2 * af::Pi * lineSpacing_Square_TiledLevel * ColumnScale / d_NoOfElements));
	//af_print(Column_ScaleFactor.slice(0))

	// Consider multiplication of complex numbers A = (a+ib); B = (c+id)	
	array realRealPart;             // ac
	array realImagPart;				// ad 
	array imagRealPart;				// bc
	array imagImagPart;				// bd

	SplitMultiplyComplex(FrFT1D_Image_3D, Column_ScaleFactor, realRealPart, realImagPart, imagRealPart, imagImagPart);
	Slice2D = sum(complex(realRealPart - imagImagPart, realImagPart + imagRealPart)); // sum(flip(FrFT_Image_X_Cube, 0) *BetaFactor);
	Slice2D_Conj = sum(complex(realRealPart + imagImagPart, imagRealPart - realImagPart));//  sum(flip(FrFT_Image_X_Cube, 0) * conjg(BetaFactor));

	/*af_print(Slice2D);
	af_print(Slice2D_Conj);*/

	Slice2D = moddims(Slice2D, d_NoOfElements, d_NoOfElements);                // Making it 2D without changing the order of data
	Slice2D_Conj = moddims(Slice2D_Conj, d_NoOfElements, d_NoOfElements);

	/*af_print(Slice2D);
	af_print(Slice2D_Conj);*/
}

void ComputeFinal_2ComplementaryLines(array& FrFTFrFT_Image2D, double ColumnScale, int d_NoOfElements, array& Line1D, array& Line1D_Conj)
{
	int N = d_NoOfElements - 1;
	array lineSpacing = array(seq(-N / 2, N / 2)).as(PRECISION_REAL);
	array lineSpacing_tiled_Y = tile(lineSpacing, 1, d_NoOfElements);
	//af_print(lineSpacing_tiled_Y);
	array lineSpacing_Square = lineSpacing_tiled_Y * lineSpacing_tiled_Y.T();
	array Column_ScaleFactor = cexp(complex(0, -2 * af::Pi * lineSpacing_Square * ColumnScale / d_NoOfElements));
	//af_print(Column_ScaleFactor)

	// Consider multiplication of complex numbers A = (a+ib); B = (c+id)	
	array realRealPart;             // ac
	array realImagPart;				// ad 
	array imagRealPart;				// bc
	array imagImagPart;				// bd

	SplitMultiplyComplex(FrFTFrFT_Image2D, Column_ScaleFactor, realRealPart, realImagPart, imagRealPart, imagImagPart);
	Line1D = sum(complex(realRealPart - imagImagPart, realImagPart + imagRealPart)); // sum(flip(FrFT_Image_X_Cube, 0) *BetaFactor);
	Line1D_Conj = sum(complex(realRealPart + imagImagPart, imagRealPart - realImagPart));//  sum(flip(FrFT_Image_X_Cube, 0) * conjg(BetaFactor));
	//af_print(Line1D);
	//af_print(Line1D_Conj);
}

array ComputeForSinglePoint3D(array& Image_3D, array& gridSpacing, int d_NoOfElements, double* desiredPoint)
{
	int N = d_NoOfElements - 1;

	// First since Columns are first dimension in ArrayFire , Y-axis scaling first
	array y_Multiplicationfactor = cexp(complex(0, -2 * Pi*desiredPoint[1] * gridSpacing / d_NoOfElements));

	array Map_y = tile(y_Multiplicationfactor, 1, d_NoOfElements, d_NoOfElements);
	//af_print(Image_3D);
	//af_print(Map_y);
	array Image_2D = moddims(sum(Image_3D*Map_y), d_NoOfElements, d_NoOfElements);
	//af_print(Image_2D);

	// X-axis 

	array x_Multiplicationfactor = cexp(-complex(0, 2 * Pi*desiredPoint[0] * gridSpacing / d_NoOfElements));

	array Map_x = tile(x_Multiplicationfactor, 1, d_NoOfElements);
	//af_print(Map_x);
	array Image_1D = sum(Image_2D*Map_x);
	//af_print(Image_1D);

	// Finally Z-axis 

	array z_Multiplicationfactor = cexp(-complex(0, 2 * Pi*desiredPoint[2] * gridSpacing / d_NoOfElements));
	Image_1D = moddims(Image_1D, d_NoOfElements, 1);
	//af_print(z_Multiplicationfactor);
	//af_print(Image_1D);
	array SinglePoint = sum(Image_1D* z_Multiplicationfactor);
	return SinglePoint;

}

void ComputeDirectBruteForce3D(array& Image_3D, int d_NoOfElements, int d_NoOfAnglesTheta, int d_NoOfAnglesPhi, array& FinalSphericalGrid)
{
	int N = d_NoOfElements - 1;
	array gridSpacing = array(seq(-N / 2, N / 2)).as(PRECISION_REAL);
	array element;
	
	double desiredPoint[3];
	

	for (int p = 0; p < d_NoOfAnglesTheta; ++p)                          // Gives the choice of Polar Slices  in the XY plane
	{
		double angleTheta = p * Pi / d_NoOfAnglesTheta;
		for (int q = 0; q < d_NoOfAnglesPhi; ++q)                       // Gives the choice of Anglular lines measured from the XY plane
		{
			double anglePhi = q * Pi / d_NoOfAnglesPhi;

			//timer::start();
			for (int n = 0; n < d_NoOfElements; ++n)                    // Gives choices of points on the angular lines
			{
				element = gridSpacing.row(n);						   // Radial distance
				double rho = element.scalar<double>();

				desiredPoint[0] = rho * cos(angleTheta)*cos(anglePhi);    // x coordinate 
				desiredPoint[1] = rho * cos(angleTheta)*sin(anglePhi);   // y coordinate 
				desiredPoint[2] = rho * sin(angleTheta);					 // z coordinate 

				//desiredPoint[0] = rho * cos(anglePhi)*cos(angleTheta);   Previous // x coordinate 
				//desiredPoint[1] = rho * cos(anglePhi)*sin(angleTheta);   // y coordinate 
				//desiredPoint[2] = rho * sin(anglePhi);					 // z coordinate 

				FinalSphericalGrid(p, q, n) = ComputeForSinglePoint3D(Image_3D, gridSpacing, d_NoOfElements, desiredPoint);
			}
			//printf("elapsed seconds: %g at iteration %d\n", timer::stop(),q);
		}
		
	}
}

array Get2DFullPolarDFT(array& Image, int d_NoOfElements, int d_NoOfAngles, int d_NoOfLevels)    // This will require the same processing as defined for 2D see the solution in 2D
{
	
	int N = d_NoOfElements - 1;
	cdouble i_cdouble = { 0, 1 };
	array i = constant(i_cdouble, 1, 1, PRECISION_COMPLEX);/* imaginary unit */

	/*-----------------------------------Preparing Padded & Tiled Images in X and Y--------------------------------------------*/
	//array Zeros = constant(0, d_NoOfElements, d_NoOfElements, PRECISION_REAL);
	array Zeros = constant(0, d_NoOfElements, d_NoOfElements, PRECISION_COMPLEX);
	array Ones = constant(1, d_NoOfElements, d_NoOfElements, PRECISION_COMPLEX);
	//af_print(Image);
	//af_print(Zeros);
	array Zero_Padded_X_Image_Transposed = join(1, Image, Zeros).T();        /* Transposing so that the rows acts as columns to be operated upon for FFTs*/
	array Zero_Padded_Y_Image = join(0, flip(Image, 0), Zeros);              /* Needed since image gathered gives +ve to -ve indexes and we require opposite for computations*/
	array Image_Tiled_X = tile(Zero_Padded_X_Image_Transposed, 1, 1, d_NoOfLevels);
	array Image_Tiled_Y = tile(Zero_Padded_Y_Image, 1, 1, d_NoOfLevels);

	/*-------------------------------------------Creating Index Cubes and Sequences----------------------------------------------------*/
	array leftSideIndexes = array(seq(0, N, 1)).as(PRECISION_REAL);
	array rightSideIndexesOnes = -1 * array(seq(1, d_NoOfElements, 1)).as(PRECISION_REAL);
	array rightSideIndexesZeros = constant(0, d_NoOfElements, 1, PRECISION_REAL);
	array rightSideIndexesN_2 = constant(N / 2, d_NoOfElements, 1, PRECISION_REAL);

	array indexedElementsEn = join(0, leftSideIndexes, flip(rightSideIndexesOnes, 0));
	array indexedElementsPre = join(0, leftSideIndexes, rightSideIndexesN_2);     /* This is for Keeping pre and post multiplication factor upper half only*/
	array indexedElementsPost = join(0, leftSideIndexes, rightSideIndexesZeros);

	array indexedElements_Tiled_En = tile(pow(indexedElementsEn, 2), 1, d_NoOfElements, d_NoOfLevels);
	array indexedElements_Tiled_PreMulti = tile(indexedElementsPre - N / 2, 1, d_NoOfElements, d_NoOfLevels);
	array indexedElements_Tiled_PostMulti = tile(indexedElementsPost, 1, d_NoOfElements, d_NoOfLevels);

	array lineSpacing = array(seq(-N / 2, N / 2)).as(PRECISION_REAL);
	array lineSpacing_tiled_Y = tile(lineSpacing, 1, d_NoOfElements);
	array lineSpacing_Square = lineSpacing_tiled_Y * lineSpacing_tiled_Y.T();
	array lineSpacing_Square_TiledLevel = tile(lineSpacing_Square, 1, 1, d_NoOfLevels);
	
	/*--------------------------Creating Alpha Cubes------------------------------------*/
	array alpha_Levels = cos(af::Pi / d_NoOfAngles * array(seq(1, d_NoOfLevels)).as(PRECISION_REAL));
	array alpha_Mods = moddims(alpha_Levels, 1, 1, d_NoOfLevels);
	array alpha_Tiled_depth = tile(alpha_Mods, 2 * d_NoOfElements, d_NoOfElements, 1);
	array alpha_Tiled = moddims(alpha_Tiled_depth, 2 * d_NoOfElements, d_NoOfElements, d_NoOfLevels);


	/*--------------------------Creating Beta Cubes---------------------------------------*/
	array beta_Levels = sin(Pi / d_NoOfAngles * array(seq(1, d_NoOfLevels)).as(PRECISION_REAL));
	array beta_Mods = moddims(beta_Levels, 1, 1, d_NoOfLevels);
	array beta_Tiled_depth = tile(beta_Mods, d_NoOfElements, d_NoOfElements, 1);
	array beta_Tiled = moddims(beta_Tiled_depth, d_NoOfElements, d_NoOfElements, d_NoOfLevels);
	

	/*-------------------Precomputing the Essential Sequence Cubes :: All complex values here --------------------*/
	array Seq_En = cexp(-tile(i, 2 * d_NoOfElements, d_NoOfElements, d_NoOfLevels) * af::Pi * indexedElements_Tiled_En * alpha_Tiled / d_NoOfElements);   /* E(n) as defined in the paper*/
	
	array subtractValues = tile(join(0, Zeros, Ones), 1, 1, d_NoOfLevels);			/* This is for  Keeping pre and post multiplication factor upper half only*/
	array PreMultiplicationFactor = cexp(tile(i, 2 * d_NoOfElements, d_NoOfElements, d_NoOfLevels) * af::Pi * indexedElements_Tiled_PreMulti * alpha_Tiled * N / d_NoOfElements) - subtractValues;
	array PostMultiplicationFactor = cexp(tile(i, 2 * d_NoOfElements, d_NoOfElements, d_NoOfLevels)* af::Pi *  indexedElements_Tiled_PostMulti * alpha_Tiled * N / d_NoOfElements) - subtractValues;
	array BetaFactor = cexp(-2 * tile(i, d_NoOfElements, d_NoOfElements, d_NoOfLevels) * af::Pi * lineSpacing_Square_TiledLevel * beta_Tiled / d_NoOfElements);

	array ZeroNinty_Factor = cexp(-2 * tile(i, d_NoOfElements, d_NoOfElements) * af::Pi * lineSpacing_Square * 1 / d_NoOfElements);

	/*--------------------Preprocessing Cubes-----------------------*/
	array Image_Tiled_X_PreMulti = Image_Tiled_X * PreMultiplicationFactor;
	array Image_Tiled_X_PreMulti_SeqEn = Image_Tiled_X_PreMulti * Seq_En;
	array Image_Tiled_Y_PreMulti = Image_Tiled_Y * PreMultiplicationFactor;
	array Image_Tiled_Y_PreMulti_SeqEn = Image_Tiled_Y_PreMulti * Seq_En;


	/*-------------------Computing Convolution--------------------*/
	int normalizationFactor = pow(2 * d_NoOfElements, 0);			// This is a scalar single value especially needed only for ArrayFire ifft
	array firstFFT_X = fft(moddims(Image_Tiled_X_PreMulti_SeqEn, 2 * d_NoOfElements, d_NoOfLevels*d_NoOfElements));
	array secondFFT_X = fft(moddims(conjg(Seq_En), 2 * d_NoOfElements, d_NoOfLevels*d_NoOfElements));
	array interim_FrFT_X = ifft(firstFFT_X * secondFFT_X) / normalizationFactor;
	array firstFFT_Y = fft(moddims(Image_Tiled_Y_PreMulti_SeqEn, 2 * d_NoOfElements, d_NoOfLevels*d_NoOfElements));
	array secondFFT_Y = secondFFT_X;
	array interim_FrFT_Y = ifft(firstFFT_Y * secondFFT_Y) / normalizationFactor;


	interim_FrFT_X = moddims(interim_FrFT_X, 2 * d_NoOfElements, d_NoOfElements, d_NoOfLevels);    /* Rearranging matrix to level based volume*/
	interim_FrFT_Y = moddims(interim_FrFT_Y, 2 * d_NoOfElements, d_NoOfElements, d_NoOfLevels);

	/*-------------------Postprocessing-----------------------------*/
	array  FrFT_Image_X = interim_FrFT_X * Seq_En * PostMultiplicationFactor;
	array  FrFT_Image_Y = interim_FrFT_Y * Seq_En * PostMultiplicationFactor;


	/*--------------------Grab only the top half elements drop overlapping------------------*/
	array FrFT_Image_X_Cube = FrFT_Image_X.rows(0, N);
	array FrFT_Image_Y_Cube = FrFT_Image_Y.rows(0, N);

	/*--------------------Finally all computations for  the  Polar Grid-----------*/
	//   Computing for all the grid expect two special indexes
	array levelSeq = array(seq(0, d_NoOfLevels - 1)).as(PRECISION_REAL);
	array finalIndexSeq1_X = 1 + levelSeq;
	array finalIndexSeq2_X = d_NoOfAngles - finalIndexSeq1_X;
	array finalIndexSeq3_Y = d_NoOfAngles / 2 - finalIndexSeq1_X;
	array finalIndexSeq4_Y = d_NoOfAngles / 2 + finalIndexSeq1_X;

	/* This transpose is Required since we will no operate in the column wise transposed axis */
	FrFT_Image_X_Cube = FrFT_Image_X_Cube.T();
	FrFT_Image_Y_Cube = FrFT_Image_Y_Cube.T();

	// Consider multiplication of complex numbers A = (a+ib); B = (c+id), This is used to reduce the computations by half 	
	array realRealPart;             // ac
	array realImagPart;				// ad
	array imagRealPart;				// bc
	array imagImagPart;				// bd

	SplitMultiplyComplex(flip(FrFT_Image_X_Cube, 0), BetaFactor, realRealPart, realImagPart, imagRealPart, imagImagPart);
	array tempSeq_X = sum(complex(realRealPart - imagImagPart, realImagPart + imagRealPart)); // sum(flip(FrFT_Image_X_Cube, 0) *BetaFactor);
	array tempSeqConj_X = sum(complex(realRealPart + imagImagPart, imagRealPart - realImagPart));//  sum(flip(FrFT_Image_X_Cube, 0) * conjg(BetaFactor));
	SplitMultiplyComplex(FrFT_Image_Y_Cube, BetaFactor, realRealPart, realImagPart, imagRealPart, imagImagPart);
	array tempSeq_Y = sum(complex(realRealPart - imagImagPart, realImagPart + imagRealPart)); //sum(FrFT_Image_Y_Cube * BetaFactor);
	array tempSeqConj_Y = sum(complex(realRealPart + imagImagPart, imagRealPart - realImagPart)); //sum(FrFT_Image_Y_Cube *conjg(BetaFactor));

	array finalSeq_X = moddims(tempSeq_X, d_NoOfElements, d_NoOfLevels).T();
	array finalSeqConj_X = moddims(tempSeqConj_X, d_NoOfElements, d_NoOfLevels).T();
	finalSeqConj_X = flip(finalSeqConj_X, 1);
	array finalSeq_Y = moddims(tempSeq_Y, d_NoOfElements, d_NoOfLevels).T();
	array finalSeqConj_Y = moddims(tempSeqConj_Y, d_NoOfElements, d_NoOfLevels).T();

	// Removing just 2 redundant computations for 45 degree case
	if (0 == remainder(d_NoOfAngles, 4))
	{
		finalIndexSeq3_Y = finalIndexSeq3_Y.rows(0, d_NoOfLevels - 2);          // Removing just the last rows from 4 structures
		finalSeq_Y = finalSeq_Y.rows(0, d_NoOfLevels - 2);
		finalIndexSeq4_Y = finalIndexSeq4_Y.rows(0, d_NoOfLevels - 2);
		finalSeqConj_Y = finalSeqConj_Y.rows(0, d_NoOfLevels - 2);
	}

	//   Computing seperately for two special indexes
	array SpecialTwoIndexes;
	array ZeroLineFrFT_Image_X_Cube;
	array NintyLineFrFT_Image_Y_Cube;
	if (PRECISION_REAL == f32)
	{
	   float zeroIndex = 0;
	   float nintyIndex = d_NoOfAngles / 2;
	    float values[] = { zeroIndex, nintyIndex };
	   SpecialTwoIndexes= array(2, 1, values);
	   ZeroLineFrFT_Image_X_Cube = FrFT_Image_Y_Cube.slice(zeroIndex).col(N / 2);
	   NintyLineFrFT_Image_Y_Cube = FrFT_Image_X_Cube.slice(zeroIndex).col(N / 2);
	}
	else
	{
		double zeroIndex = 0;
		double nintyIndex = d_NoOfAngles / 2;
		double values[] = { zeroIndex, nintyIndex };
		SpecialTwoIndexes = array(2, 1, values);
		ZeroLineFrFT_Image_X_Cube = FrFT_Image_Y_Cube.slice(zeroIndex).col(N / 2);
		NintyLineFrFT_Image_Y_Cube = FrFT_Image_X_Cube.slice(zeroIndex).col(N / 2);
	}

	

	array DFTZeroLine = sum(tile(ZeroLineFrFT_Image_X_Cube, 1, d_NoOfElements) *ZeroNinty_Factor);
	array DFTNinetyLine = sum(tile(flip(NintyLineFrFT_Image_Y_Cube, 0), 1, d_NoOfElements) *ZeroNinty_Factor);
	array SpecialTwoLines = join(0, DFTZeroLine, DFTNinetyLine);


	array UnsortedIndexes = join(0, join(0, join(0, join(0, finalIndexSeq1_X, finalIndexSeq2_X), finalIndexSeq3_Y), finalIndexSeq4_Y), SpecialTwoIndexes);
	array tiledUnsortedIndexes = tile(UnsortedIndexes, 1, d_NoOfElements);
	array UnsortedPolarGrid = join(0, join(0, join(0, join(0, finalSeq_X, finalSeqConj_X), finalSeq_Y), finalSeqConj_Y), SpecialTwoLines);


	array FinalPolarGridReal;// = constant(0, d_NoOfElements, d_NoOfAngles, PRECISION_COMPLEX);
	array Output_Keys_Sorted;
	sort(Output_Keys_Sorted, FinalPolarGridReal, tiledUnsortedIndexes, real(UnsortedPolarGrid));

	array FinalPolarGridImag;// = constant(0, d_NoOfElements, d_NoOfAngles, PRECISION_COMPLEX);
	array Output_Keys_Sorted2;
	sort(Output_Keys_Sorted2, FinalPolarGridImag, tiledUnsortedIndexes, imag(UnsortedPolarGrid));


	return complex(FinalPolarGridReal, FinalPolarGridImag);
}

static void ComputeFastSphericalPolarFourierTransform(array& Image3D, array& FinalSphericalGrid, int d_NoOfElements, int d_NoOfAnglesTheta, int d_NoOfLevelsTheta, int d_NoOfAnglesPhi, int d_NoOfLevelsPhi)
{
	array ReorderedImage_OperateColumns;                 // Reordered 3D Image at the start of the FrFT operation  
	array FrFT1D_Uniform_Image3D;                        // 3D Image obtained after passing 1st stage of FrFT, all columns have now been operated on
	array Silce2D, Slice2D_Conj;                         // 2D Slices obtained at the 2nd stage of FrFT
	array Line1, Line1_Conj, Line2, Line2_Conj;			 // 1D Lines at the end of the 3rd stage of FrFT, which are final	
	double alpha_factor, beta_factor, gamma_factor;      // Scaling factors in X-axis , Y-axis and Z-axis respectively, which change depending on the computation of the block 
	
	int N = d_NoOfElements - 1;
	GlobalArraysComputeInitialize(d_NoOfElements);
	for (int q = 1; q <= d_NoOfLevelsPhi; ++q)                // Splitting as defined in the paper (See the schematic of spherical divisions)
	{
		double anglePhi = q * Pi / d_NoOfAnglesPhi;

		//timer::start();
		for (int p = 1; p <= d_NoOfLevelsTheta; ++p)       // X oriented pair of Polar slices, ( q, d_NoOfAnglesPhi - q)
		{
			double angleTheta = p  * Pi / d_NoOfAnglesTheta;

			/**************************************************************************************************************************************************************************/
			// XX block  -- Concentric rectangles in YZ tiled along X -axis
			alpha_factor = cos(angleTheta) * cos(anglePhi);											         // Scaling needed as defined in the paper
			beta_factor = cos(angleTheta) * sin(anglePhi);
			gamma_factor = sin(angleTheta);

			ReorderedImage_OperateColumns = reorder(Image3D, X_AXIS, Y_AXIS, Z_AXIS);                                //  1st Swap, X-axis as columns now
			//af_print(Image3D.T().slices(0, 2));
			//af_print(ReorderedImage_OperateColumns.slices(0,2));
			//ReorderedImage_OperateColumns = reorder(ReorderedImage_OperateColumns, X_AXIS, Y_AXIS, Z_AXIS );            // Restoring the image as it was
			//af_print(ReorderedImage_OperateColumns.slices(0,2));

			FrFT1D_Uniform_Image3D = ComputeNonVectorColumnwise_FrFTUniform(ReorderedImage_OperateColumns, alpha_factor, d_NoOfElements); // This is the most expensive operation of order (N+1)^3 log(N+1)
			//af_print(FrFT1D_Uniform_Image3D)
			if (q == 1 && p == 1)         // Computing Polar Slice at anglephi = 90 , only need to be computed once
			{
				array Central2D_YZSlice2D = FrFT1D_Uniform_Image3D(N / 2, span, span);
				Central2D_YZSlice2D = moddims(Central2D_YZSlice2D, d_NoOfElements, d_NoOfElements);
				//af_print(Central2D_YZSlice2D);

				array properOrientedSlice = flip(Central2D_YZSlice2D.T(), 0);          // Verified match 
				//af_print(properOrientedSlice);
				array Polar2D = Get2DFullPolarDFT(properOrientedSlice, d_NoOfElements, d_NoOfAnglesTheta, d_NoOfLevelsTheta);
				//af_print(FinalSphericalGrid(span, d_NoOfAnglesPhi / 2, span));
				Polar2D = moddims(Polar2D, d_NoOfAnglesTheta, 1, d_NoOfElements);
				//af_print(Polar2D);
				FinalSphericalGrid(span, d_NoOfAnglesPhi / 2, span) = Polar2D;
				//af_print(BruteForceSphericalGrid(span, d_NoOfAnglesPhi / 2, span) - FinalSphericalGrid(span, d_NoOfAnglesPhi / 2, span));
			}

			//af_print(FrFT1D_Uniform_Image3D.slice(0))
			//af_print((reorder(FrFT1D_Uniform_Image3D, Z_AXIS, Y_AXIS, X_AXIS)).slice(0))
			//af_print(FrFT1D_Uniform_Image3D);
			//af_print(reorder(FrFT1D_Uniform_Image3D, Y_AXIS, X_AXIS, Z_AXIS));
			ComputeNonVector_NonUniformFrFT_SingleColumn(reorder(FrFT1D_Uniform_Image3D, X_AXIS, Y_AXIS, Z_AXIS), beta_factor, d_NoOfElements, Silce2D, Slice2D_Conj); //  2nd Swap, Y-axis as columns
			//af_print(Silce2D);
			//af_print(Slice2D_Conj);

			// Gather 4 lines
			ComputeFinal_2ComplementaryLines(Silce2D.T(), gamma_factor, d_NoOfElements, Line1, Line1_Conj);          // 3rd Swap , Z axis as columns
			ComputeFinal_2ComplementaryLines(Slice2D_Conj.T(), gamma_factor, d_NoOfElements, Line2, Line2_Conj);      // 3rd Swap , Z axis as columns

			//af_print(Line1);
			//af_print(Line1_Conj); 
			//af_print(Line2);
			//af_print(Line2_Conj);

			FinalSphericalGrid(p, q, span) = Line1;
			FinalSphericalGrid(p, d_NoOfAnglesPhi - q, span) = conjg(Line2_Conj);
			FinalSphericalGrid(d_NoOfAnglesTheta - p, q, span) = conjg(Line1_Conj);											// 1 More Swap to match !
			FinalSphericalGrid(d_NoOfAnglesTheta - p, d_NoOfAnglesPhi - q, span) = Line2;

			/*printf("\n My Special Solution 1st XX block\n");
			af_print(moddims(FinalSphericalGrid(p, q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(p, d_NoOfAnglesPhi - q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta - p, q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta - p, d_NoOfAnglesPhi - q, span), 1, d_NoOfElements));
			printf("\n");*/
			/**************************************************************************************************************************************************************************/

			if (p * 180 / d_NoOfAnglesTheta == 45)       // At 45 degrees its redundant to compute in the next block so skip it
				continue;


			/**************************************************************************************************************************************************************************/
			// XZ block  -- Concentric rectangles in XY tiled along Z -axis
			alpha_factor = sin(angleTheta) * cos(anglePhi);											         // Scaling needed as defined in the paper
			beta_factor = sin(angleTheta) * sin(anglePhi);
			gamma_factor = cos(angleTheta);

			ReorderedImage_OperateColumns = reorder(Image3D, Z_AXIS, Y_AXIS, X_AXIS);                                   //  1st Swap, Z-axis elements in columns now
			//af_print(ReorderedImage_OperateColumns.slices(0,2));
			//ReorderedImage_OperateColumns = reorder(ReorderedImage_OperateColumns, X_AXIS, Z_AXIS, Y_AXIS );           // Restoring the image as it was
			//af_print(ReorderedImage_OperateColumns.slices(0, 2));
			FrFT1D_Uniform_Image3D = ComputeNonVectorColumnwise_FrFTUniform(ReorderedImage_OperateColumns, gamma_factor, d_NoOfElements); // This is the most expensive operation of order (N+1)^3 log(N+1)
			//af_print(FrFT1D_Uniform_Image3D);

			if (q == 1 && p == 1)         // Computing One special Polar Slice at angletheta = 0 , only need to be computed once
			{
				array Central2D_XYSlice2D = FrFT1D_Uniform_Image3D(N / 2, span, span);
				Central2D_XYSlice2D = moddims(Central2D_XYSlice2D, d_NoOfElements, d_NoOfElements);
				//af_print(Central2D_XYSlice2D);

				array properOrientedSlice = flip(Central2D_XYSlice2D, 0);          // Verified match 
				array Polar2D = Get2DFullPolarDFT(properOrientedSlice, d_NoOfElements, d_NoOfAnglesPhi, d_NoOfLevelsPhi);
				Polar2D = moddims(Polar2D, 1, d_NoOfAnglesPhi, d_NoOfElements);
				//af_print(Polar2D);
				FinalSphericalGrid(0, span, span) = Polar2D;
				//af_print(BruteForceSphericalGrid(0, span, span) - FinalSphericalGrid(0, span, span));
			}


			ComputeNonVector_NonUniformFrFT_SingleColumn(FrFT1D_Uniform_Image3D.T(), beta_factor, d_NoOfElements, Silce2D, Slice2D_Conj); //  2nd Swap, Y-axis as columns

			// Gather 4 lines
			ComputeFinal_2ComplementaryLines(Silce2D.T(), alpha_factor, d_NoOfElements, Line1, Line1_Conj);          // 3rd Swap , X axis as columns
			ComputeFinal_2ComplementaryLines(Slice2D_Conj.T(), alpha_factor, d_NoOfElements, Line2, Line2_Conj);     // Swap Lines ! // 3rd Swap , Z axis as columns

			FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, q, span) = Line1;
			FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, d_NoOfAnglesPhi - q, span) = (Line1_Conj);

			FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, q, span) = (Line2_Conj);
			FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, d_NoOfAnglesPhi - q, span) = (Line2);

			/*printf("\n My Special Solution 1st XZ block\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, d_NoOfAnglesPhi - q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, d_NoOfAnglesPhi - q, span), 1, d_NoOfElements));
			printf("\n");*/

			/**************************************************************************************************************************************************************************/

			if (q * 180 / d_NoOfAnglesPhi == 45)       // At 45 degrees its redundant to compute same in the next block  so skip it
				continue;

			/**************************************************************************************************************************************************************************/
			// YZ block                                                   // Y oriented pair of Polar slices, (d_NoOfAnglesPhi / 2 - q, d_NoOfAnglesPhi / 2 + q)
			alpha_factor = sin(angleTheta) * sin(anglePhi);				 // Scaling needed as defined in the paper
			beta_factor = sin(angleTheta) * cos(anglePhi);
			gamma_factor = cos(angleTheta);

			// Reusing the previously computed FrFT 1D Uniform scaling since the gamma_factor for both XZ and YZ block is the same

			//ReorderedImage_OperateColumns = reorder(Image3D, Z_AXIS, Y_AXIS, X_AXIS);                                   //  1st Swap, Z-axis elements in columns now
			//FrFT1D_Uniform_Image3D = ComputeNonVectorColumnwise_FrFTUniform(ReorderedImage_OperateColumns, gamma_factor, d_NoOfElements); // This is the most expensive operation of order (N+1)/4*(N+1)^3 log(N+1)
			ComputeNonVector_NonUniformFrFT_SingleColumn(FrFT1D_Uniform_Image3D.T(), beta_factor, d_NoOfElements, Silce2D, Slice2D_Conj); //  2nd Swap, Y-axis as columns

			// Gather 4 lines
			ComputeFinal_2ComplementaryLines(Silce2D.T(), alpha_factor, d_NoOfElements, Line1, Line1_Conj);          // 3rd Swap , X axis as columns
			ComputeFinal_2ComplementaryLines(Slice2D_Conj.T(), alpha_factor, d_NoOfElements, Line2, Line2_Conj);     // Swap Lines ! // 3rd Swap , Z axis as columns

			FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, d_NoOfAnglesPhi / 2 - q, span) = Line1;
			FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, d_NoOfAnglesPhi / 2 + q, span) = Line1_Conj;

			FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, d_NoOfAnglesPhi / 2 - q, span) = Line2_Conj;
			FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, d_NoOfAnglesPhi / 2 + q, span) = Line2;

			/*printf("\n My Special Solution 1st YZ block\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, d_NoOfAnglesPhi / 2 - q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 - p, d_NoOfAnglesPhi / 2 + q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, d_NoOfAnglesPhi / 2 - q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta / 2 + p, d_NoOfAnglesPhi / 2 + q, span), 1, d_NoOfElements));
			printf("\n");*/
			/**************************************************************************************************************************************************************************/

		}
		//printf("elapsed seconds: %g at iteration %d X-oriented Polar Slices\n", timer::stop(),q);

		if (q * 180 / d_NoOfAnglesPhi == 45)       // At 45 degrees its redundant to compute same in the next block  so skip it
			continue;
		//timer::start();

		for (int p = 1; p <= d_NoOfLevelsTheta; ++p)        // Y oriented pair of Polar slices, (d_NoOfAnglesPhi / 2 - q, d_NoOfAnglesPhi / 2 + q)
		{
			double angleTheta = p  * Pi / d_NoOfAnglesTheta;

			/**************************************************************************************************************************************************************************/
			// YY block
			alpha_factor = cos(angleTheta) * sin(anglePhi);				 // Scaling needed as defined in the paper
			beta_factor = cos(angleTheta) * cos(anglePhi);
			gamma_factor = sin(angleTheta);

			ReorderedImage_OperateColumns = Image3D;                                         // No swap directly we can operate on columns
			FrFT1D_Uniform_Image3D = ComputeNonVectorColumnwise_FrFTUniform(ReorderedImage_OperateColumns, beta_factor, d_NoOfElements); // This is the most expensive operation of order (N+1)^3 log(N+1)

			if (q == 1 && p == 1)         // Computing Polar Slice at anglephi = 0, only need to be computed once
			{ 
				array Central2D_YZSlice2D = FrFT1D_Uniform_Image3D(N / 2, span, span);
				Central2D_YZSlice2D = moddims(Central2D_YZSlice2D, d_NoOfElements, d_NoOfElements);
				//af_print(Central2D_YZSlice2D);

				array properOrientedSlice = flip(Central2D_YZSlice2D.T(), 0);          // Verified match  
				array Polar2D = Get2DFullPolarDFT(properOrientedSlice, d_NoOfElements, d_NoOfAnglesTheta, d_NoOfLevelsTheta);
				Polar2D = moddims(Polar2D, d_NoOfAnglesTheta, 1, d_NoOfElements);
				//af_print(Polar2D);
				FinalSphericalGrid(span, 0, span) = Polar2D;
				//af_print(BruteForceSphericalGrid(span, 0, span) - FinalSphericalGrid(span, 0, span));


				// Special Operation
				//af_print(Polar2D(d_NoOfAnglesTheta / 2, span));
				array SpecialLineZ = Polar2D(d_NoOfAnglesTheta / 2, span);
				//af_print(SpecialLineZ);
				array TiledZ = tile(SpecialLineZ, 1, d_NoOfAnglesPhi, 1);
				//af_print(TiledZ);
				//af_print(FinalSphericalGrid(d_NoOfAnglesTheta / 2, span, span));
				FinalSphericalGrid(d_NoOfAnglesTheta / 2, span, span) = TiledZ;        // This line is common to all grids !!! VERY VERY Special
			}

			ComputeNonVector_NonUniformFrFT_SingleColumn(FrFT1D_Uniform_Image3D.T(), alpha_factor, d_NoOfElements, Silce2D, Slice2D_Conj); //  2nd Swap, X-axis as columns now

			// Gather 4 lines
			ComputeFinal_2ComplementaryLines(Silce2D.T(), gamma_factor, d_NoOfElements, Line1, Line1_Conj);          // 3rd Swap , Z axis as columns
			ComputeFinal_2ComplementaryLines(Slice2D_Conj.T(), gamma_factor, d_NoOfElements, Line2, Line2_Conj);     // Swap Lines ! // 3rd Swap , Z axis as columns


			FinalSphericalGrid(p, d_NoOfAnglesPhi / 2 - q, span) = Line1;
			FinalSphericalGrid(p, d_NoOfAnglesPhi / 2 + q, span) = Line2;

			FinalSphericalGrid(d_NoOfAnglesTheta - p, d_NoOfAnglesPhi / 2 - q, span) = conjg(Line1_Conj);
			FinalSphericalGrid(d_NoOfAnglesTheta - p, d_NoOfAnglesPhi / 2 + q, span) = conjg(Line2_Conj);

			/*printf("\n My Special Solution 1st YY block\n");
			af_print(moddims(FinalSphericalGrid(p, d_NoOfAnglesPhi / 2 - q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(p, d_NoOfAnglesPhi / 2 + q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta - p, d_NoOfAnglesPhi / 2 - q, span), 1, d_NoOfElements));
			printf("\n");
			af_print(moddims(FinalSphericalGrid(d_NoOfAnglesTheta - p, d_NoOfAnglesPhi / 2 + q, span), 1, d_NoOfElements));
			printf("\n");*/
			/**************************************************************************************************************************************************************************/

			//if (p * 180 / d_NoOfAnglesTheta == 45)       // At 45 degrees its redundant to compute in the next block       so skip it
			//	continue;

			
		}
		//printf("elapsed seconds: %g at iteration %d Y-oriented Polar Slices\n", timer::stop(),q);
	}
	
}

static int d_NoOfElements = 9;		/* IT is always odd */

static void bench()
{
	array Image3D = randu(d_NoOfElements, d_NoOfElements, d_NoOfElements, PRECISION_REAL);
	int d_NoOfAnglesTheta = d_NoOfElements + 13;       /*Always even*/
	int d_NoOfAnglesPhi = d_NoOfElements + 29;        /*Always even*/
	int d_NoOfLevelsTheta = ceil(float(d_NoOfAnglesTheta - 2) / 4);          // This is where partitioning begins
	int d_NoOfLevelsPhi = ceil(float(d_NoOfAnglesPhi - 2) / 4);

	array FinalSphericalGrid = constant(0, d_NoOfAnglesTheta, d_NoOfAnglesPhi, d_NoOfElements, PRECISION_COMPLEX);      // Angle Phi vs. Polar slices : Each Polar slice has No of elevation angle Theta vs. Radial data
	ComputeFastSphericalPolarFourierTransform(Image3D, FinalSphericalGrid, d_NoOfElements, d_NoOfAnglesTheta, d_NoOfLevelsTheta, d_NoOfAnglesPhi, d_NoOfLevelsPhi);

}

static void MakeFileEntriesAsInputs()
{
	std::ofstream  myfile;
	myfile.open("Input.csv");
	//myfile << "InputArguments\n";
	for (int noOfElements = 15; noOfElements < 168; noOfElements = noOfElements + 2)
	{
		myfile << "{Back}{Back}{Back}" << noOfElements <<  "\n";
	}
	myfile.close();
}

static void MakeFileEntriesAsResults(int noOfElements)
{
	std::ofstream  myfile;
	myfile.open("Results.csv", std::ios::app);
	//myfile << "Timing with microsecond accurate clock (GPU) in seconds\n \n \n ";
	//myfile << "(N+1), times(seconds)\n";
	
	
	//for (int noOfElements = 15; noOfElements < 178; noOfElements = noOfElements + 2)
	{
		af::deviceGC();
		d_NoOfElements = noOfElements;
		printf("%dx%dx%d random number generated 3D Image and its Spherical Polar FFT \n", d_NoOfElements, d_NoOfElements, d_NoOfElements);
		double time_s = timeit(bench); // seconds
		printf("Timing with microsecond accurate clock (GPU) %f seconds\n", time_s);
		myfile << d_NoOfElements << "," << time_s << "\n";
	}

	myfile.close();
}



BOOL KillProcessByName(char *szProcessToKill){
	HANDLE hProcessSnap;
	HANDLE hProcess;
	PROCESSENTRY32 pe32;
	DWORD dwPriorityClass;

	hProcessSnap = CreateToolhelp32Snapshot(TH32CS_SNAPPROCESS, 0);  // Takes a snapshot of all the processes

	if (hProcessSnap == INVALID_HANDLE_VALUE){
		return(FALSE);
	}

	pe32.dwSize = sizeof(PROCESSENTRY32);

	if (!Process32First(hProcessSnap, &pe32)){
		CloseHandle(hProcessSnap);
		return(FALSE);
	}

	do{
		if (!strcmp(pe32.szExeFile, szProcessToKill)){    //  checks if process at current position has the name of to be killed app
			hProcess = OpenProcess(PROCESS_TERMINATE, 0, pe32.th32ProcessID);  // gets handle to process
			TerminateProcess(hProcess, 0);   // Terminate process by handle
			CloseHandle(hProcess);  // close the handle
		}
	} while (Process32Next(hProcessSnap, &pe32));  // gets next member of snapshot

	CloseHandle(hProcessSnap);  // closes the snapshot handle
	return(TRUE);
}


 
int main(int argc, char* argv[]) 
{
	try {
		  
		info();
		//int total_device = devicecount();
		//printf("This computer has device count %d \n", total_device);

		//MakeFileEntriesAsInputs();
		//int noOfElements;
		//sscanf(argv[1], "%d", &noOfElements);
		//MakeFileEntriesAsResults(noOfElements);
		
		//KillProcessByName("cmd.exe");
		
		/******************************* Check the result ! *************************************************/
		/***************************************************************************************************/
		/***************************************************************************************************/
		/***************************************************************************************************/
		/* With d_NoOfElements = 129 which seems moderate size of 3D image, this is the result of computations: 
		ArrayFire v3.0.0 (CUDA, 64 - bit Windows, build 86426db17
		Platform: CUDA Toolkit 7, Driver : 0.00
		[0] GeForce GTX 760, 2048 MB, CUDA Compute 3.0
		This computer has device count 1
		129x129x129 random number generated 3D Image and its Spherical Polar FFT --- These timings are only for CPU sides, accurate timing with GPU clock is lower 
		Total time with brute force O((N + 1) ^ 6), elapsed seconds : 14667.2                                     --- This is approx 4 hours !!!
		Total time for my Special solution O(3 x K/4 x M/4 x (N+1)^3/2 x log2(N+1)), elapsed seconds : 410        --- This is approx 6.83 mins !!!!!!!!!!!!

		/*******************************************************************************************************
		with GeForce GTX Titan X we have ...
		Total time for my Special solution O(3 x K/4 x M/4 x (N+1)^3/2 x log2(N+1)), elapsed seconds : 254       --- This is approx 4.233 mins !!!!!!!!!!!!

		For d_NoOfElements = 133, elapsed seconds = 419 seconds
		Error between the direct brute force computation and the fast exact solution
		abs(sum(sum(sum(BruteForceSphericalGrid - FinalSphericalGrid))))[1 1 1 1]
		0.0
		hit[enter]...   */
		/***************************************************************************************************/
		/***************************************************************************************************/
		/***************************************************************************************************/

		 
		/******************************Create DATA********************************************/
		 
		int total_device = devicecount();
		printf("This computer has device count %d \n", total_device);
		printf("%dx%dx%d random number generated 3D Image and its Spherical Polar FFT \n", d_NoOfElements, d_NoOfElements, d_NoOfElements);
		 
		array Image3D(d_NoOfElements, d_NoOfElements, d_NoOfElements, ThreeDImageValues);   // Gather image from the stored values  for testing purposes
//		array Image3D = randu(d_NoOfElements, d_NoOfElements, d_NoOfElements,  PRECISION_REAL);
		Image3D = Image3D.T();
		af_print(Image3D(span,span,0));

		/*****************Learning to display 3 axis data*****************************/
		printf("This is Y -axis !\n");
		af_print(Image3D.slice(0).col(0));			// Y-Axis
		printf("This is X -axis !\n");
		af_print(Image3D.slice(0).row(0));			// X-axis
		printf("This is Z -axis !\n");
		af_print(Image3D(0, 0, span));                // Z-axis
		printf("Successfully completed the computations of 3D DFT on a full Polar Grid !\n");

		/*********************************Processing FFTs along a given axis ************************************************/
		array FFT_Y = fft(Image3D);
		printf("This is Y -axis FFTs !\n");
		af_print(fft(Image3D.slice(0).col(0)));			       // Y-axis
		af_print(FFT_Y.slice(0).col(0))

		array ReorderedImage = reorder(Image3D, X_AXIS, Y_AXIS, Z_AXIS);
		af_print(ReorderedImage.slices(1, 3));
		array TransposedImage3d = Image3D.T();
		af_print(TransposedImage3d.slices(1, 3));
		array FFT_X = fft(TransposedImage3d);
		printf("This is X -axis FFTs !\n");
		af_print(fft(TransposedImage3d.slice(0).col(0)));		// X-axis
		af_print(FFT_X.slice(0).col(0))

		array ZAxis_TransposedImage3D = reorder(Image3D, 2, 0, 1) ;
		array FFT_Z = fft(ZAxis_TransposedImage3D);
		printf("This is Z -axis FFTs !\n");
		af_print(fft(ZAxis_TransposedImage3D.slice(0).col(0)));		// X-axis
		af_print(FFT_Z.slice(0).col(0))
		 
		int N = d_NoOfElements - 1; /* it is always even as described in the paper*/
		int d_NoOfAnglesTheta = d_NoOfElements + 1;       /*Always even*/
		int d_NoOfAnglesPhi = d_NoOfElements + 5;        /*Always even*/

		
		/***************************          Brute Force solution  order   (N+1)^6       *************************************************/
		 //start timer
		timer::start();
		array BruteForceSphericalGrid = constant(0, d_NoOfAnglesTheta, d_NoOfAnglesPhi, d_NoOfElements, PRECISION_COMPLEX);      // Angle Phi vs. Polar slices : Each Polar slice has No of angle Theta vs. Radial data
		ComputeDirectBruteForce3D(Image3D, d_NoOfElements, d_NoOfAnglesTheta, d_NoOfAnglesPhi, BruteForceSphericalGrid); // Filling up the solution data
		printf("Total time with brute force O((N+1)^6), elapsed seconds: %g\n", timer::stop());
		af_print(BruteForceSphericalGrid);

		int p__ = 1;																				  // Angle theta = 18, Phi = 18
		int q__ = 1;
		printf("\n Brute Force Solution 1st XX block\n");
		af_print(moddims(BruteForceSphericalGrid(p__, q__, span), 1, d_NoOfElements));										   // Line 1
		printf("\n");
		af_print(moddims(BruteForceSphericalGrid(p__, d_NoOfAnglesPhi - q__, span), 1, d_NoOfElements));					   // Line 2
		printf("\n");
		af_print(moddims(BruteForceSphericalGrid(d_NoOfAnglesTheta - p__, q__, span), 1, d_NoOfElements));					   // Line 3
		printf("\n");
		af_print(moddims(BruteForceSphericalGrid(d_NoOfAnglesTheta - p__, d_NoOfAnglesPhi - q__, span), 1, d_NoOfElements));   // Line 4
		printf("\n");

		printf("\n Brute Force Solution 1st XZ block\n");
		af_print(moddims(BruteForceSphericalGrid(d_NoOfAnglesTheta / 2 - p__, q__, span), 1, d_NoOfElements));					    // Line 1
		printf("\n");
		af_print(moddims(BruteForceSphericalGrid(d_NoOfAnglesTheta / 2 - p__, d_NoOfAnglesPhi - q__, span), 1, d_NoOfElements));	// Line 2
		printf("\n");
		af_print(moddims(BruteForceSphericalGrid(d_NoOfAnglesTheta / 2 + p__, q__, span), 1, d_NoOfElements));						// Line 3
		printf("\n");
		af_print(moddims(BruteForceSphericalGrid(d_NoOfAnglesTheta / 2 + p__, d_NoOfAnglesPhi - q__, span), 1, d_NoOfElements));    // Line 4
		printf("\n");

		
		printf("\n Brute Force Solution 1st YY block\n");
		af_print(moddims(BruteForceSphericalGrid(p__, d_NoOfAnglesPhi / 2 - q__, span), 1, d_NoOfElements));					   // Line 1
		printf("\n");
		af_print(moddims(BruteForceSphericalGrid(p__, d_NoOfAnglesPhi / 2 + q__, span), 1, d_NoOfElements));					   // Line 2
		printf("\n");
		af_print(moddims(BruteForceSphericalGrid(d_NoOfAnglesTheta - p__, d_NoOfAnglesPhi / 2 - q__, span), 1, d_NoOfElements));   // Line 3
		printf("\n");
		af_print(moddims(BruteForceSphericalGrid(d_NoOfAnglesTheta - p__, d_NoOfAnglesPhi / 2 + q__, span), 1, d_NoOfElements));   // Line 4
		printf("\n");
		

		printf("\n Brute Force Solution 1st YZ block\n");
		af_print(moddims(BruteForceSphericalGrid(d_NoOfAnglesTheta / 2 - p__, d_NoOfAnglesPhi / 2 - q__, span), 1, d_NoOfElements));	// Line 1
		printf("\n");
		af_print(moddims(BruteForceSphericalGrid(d_NoOfAnglesTheta / 2 - p__, d_NoOfAnglesPhi / 2 + q__, span), 1, d_NoOfElements));	// Line 2
		printf("\n");
		af_print(moddims(BruteForceSphericalGrid(d_NoOfAnglesTheta / 2 + p__, d_NoOfAnglesPhi / 2 - q__, span), 1, d_NoOfElements));     // Line 3
		printf("\n");
		af_print(moddims(BruteForceSphericalGrid(d_NoOfAnglesTheta / 2 + p__, d_NoOfAnglesPhi / 2 + q__, span), 1, d_NoOfElements));	 // Line 4
		printf("\n");

		/******************************Begin Actual Spherical DFT****************************/
		double time_s = timeit(bench); // seconds
		printf("Timing with microsecond accurate clock (GPU) %f seconds\n", time_s);
		 
		timer::start();
		int d_NoOfLevelsTheta = ceil(float(d_NoOfAnglesTheta - 2) / 4);          // This is where partitioning begins
		int d_NoOfLevelsPhi = ceil(float(d_NoOfAnglesPhi - 2) / 4);

		array FinalSphericalGrid = constant(0, d_NoOfAnglesTheta, d_NoOfAnglesPhi, d_NoOfElements, PRECISION_COMPLEX);      // Angle Phi vs. Polar slices : Each Polar slice has No of elevation angle Theta vs. Radial data
		ComputeFastSphericalPolarFourierTransform(Image3D, FinalSphericalGrid, d_NoOfElements, d_NoOfAnglesTheta, d_NoOfLevelsTheta, d_NoOfAnglesPhi, d_NoOfLevelsPhi);
		printf("Total time for my Special solution O( 3 x K/4 x M/4 x (N+1)^3/2 x log2(N+1) ), elapsed seconds: %g\n", timer::stop());
		
		
		printf("\nBrute Force Special Grid\n");
		af_print(BruteForceSphericalGrid);
		printf("\nMy Special Solution\n");
		af_print(FinalSphericalGrid);

		
		af_print(FinalSphericalGrid(span, 0, span) - BruteForceSphericalGrid(span, 0, span));            // This is YZ slice anglePhi == 0
		
		af_print(FinalSphericalGrid(span, d_NoOfAnglesPhi / 2, span) - BruteForceSphericalGrid(span, d_NoOfAnglesPhi / 2, span));            // This is YZ slice anglePhi == 90

		af_print(FinalSphericalGrid(0, span, span) - BruteForceSphericalGrid(0, span, span));                  // This is the XY slice when angleTheta == 0 

		af_print(FinalSphericalGrid(d_NoOfAnglesTheta / 2, span, span) - BruteForceSphericalGrid(d_NoOfAnglesTheta / 2, span, span));
		

		printf("\nError between the direct brute force computation and the fast exact solution \n");
		af_print(BruteForceSphericalGrid - FinalSphericalGrid); 
		af_print(abs(sum(sum(sum(BruteForceSphericalGrid - FinalSphericalGrid)))));

	}
	catch (af::exception& e) {
		fprintf(stderr, "%s\n", e.what());
		throw;
	}

//#ifdef WIN32 // pause in Windows
//	if (!(argc == 2 && argv[1][0] == '-')) {
//		printf("hit [enter]...");
//		getchar();
//	}
//#endif
	return 0;
	
	
}





