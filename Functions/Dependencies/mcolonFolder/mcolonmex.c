/**************************************************************************
 * MATLAB mex function mcolonmex
 *
 * R = mcolonmex(start,step,length)
 *
 * INPUTS: start, step, length are arrays of the same length
 *         start and step must have the same class.
 *         length is double array
 * R : concatenate of colons [I1, I2, ... ] where: 
 *    Ii :=[ start(i)+step(i)*j ] for j=1, ..., length(i).
 *
 * Compile:
 *      mex mcolonmex.c on 32-bit
 *      mex -largeArrayDims mcolonmex.c on 64-bit
 *
 * MEX engine used by MCOLON.M
 *
 * Author: Bruno Luong <brunoluong@yahoo.com>
 * Date: 31-Dec-2010
 *
**************************************************************************/

#include "mex.h"

/* Define correct type depending on platform
  You might have to modify here depending on your compiler */
#if defined(_MSC_VER) || defined(__BORLANDC__)
typedef __int64 int64;
typedef __int32 int32;
typedef __int16 int16;
typedef __int8 int08;
typedef unsigned __int64 uint64;
typedef unsigned __int32 uint32;
typedef unsigned __int16 uint16;
typedef unsigned __int8 uint08;
#else /* LINUX + LCC, CAUTION: not tested by the author */
typedef long long int int64;
typedef long int int32;
typedef short int16;
typedef char int08;
typedef unsigned long long int uint64;
typedef unsigned long int uint32;
typedef unsigned short uint16;
typedef unsigned char uint08;
#endif
        
#define FILENGINE(pA, pS, pR, a, s, Type) \
    pA = (Type*)mxGetData(A); \
    pS = (Type*)mxGetData(STEP); \
    pL = PrL; \
    lgt = 0; \
    for (i=0; i<n; i++) lgt += (mwSize)(*(pL++)); \
    R = mxCreateNumericMatrix(1, 0, ClassID, mxREAL); \
    mxSetN(R, lgt); \
    pR = (Type*)mxMalloc(sizeof(Type)*lgt); \
    mxSetData(R, (void*)pR); \
    pL = PrL; \
    for (i=0; i<n; i++) { \
        lgt = (mwSize)(*(pL++)); \
        a = *(pA++); \
        s = *(pS++); \
        if (s == 1) \
            for (j=0; j<lgt; j++) *(pR++) = a + (Type)j; \
        else \
            for (j=0; j<lgt; j++) *(pR++) = a + s*(Type)j; \
    }
        
#define A prhs[0]
#define STEP prhs[1]
#define LENGTH prhs[2]
#define R plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    mxClassID ClassID;
    mwSize i, j, m, n, lgt;
    
    double *PrA, *PrS, *PrL, *pL;
    
    double *pAdouble, *pSdouble, *pRdouble, adouble, sdouble;
    float *pAsingle, *pSsingle, *pRsingle, asingle, ssingle;
    int64 *pAi64, *pSi64, *pRi64, ai64, si64;
    int32 *pAi32, *pSi32, *pRi32, ai32, si32;
    int16 *pAi16, *pSi16, *pRi16, ai16, si16;
    int08 *pAi08, *pSi08, *pRi08, ai08, si08;
    uint64 *pAu64, *pSu64, *pRu64, au64, su64;
    uint32 *pAu32, *pSu32, *pRu32, au32, su32;
    uint16 *pAu16, *pSu16, *pRu16, au16, su16;
    uint08 *pAu08, *pSu08, *pRu08, au08, su08;
    
    if( nrhs != 3 ) mexErrMsgTxt("Three arguments required.");
    
    m = mxGetM(A);
    n = mxGetN(A);
    n *= m;
    
    /* Get data pointers */
    PrL = mxGetPr(LENGTH);
    
    /* Get class of input matrix */
    ClassID = mxGetClassID(A);
    
    switch (ClassID) {
        case mxDOUBLE_CLASS:
            FILENGINE(pAdouble, pSdouble, pRdouble, adouble, sdouble, double);
            return;
        case mxSINGLE_CLASS:
            FILENGINE(pAsingle, pSsingle, pRsingle, asingle, ssingle, float);
            return;
        case mxINT64_CLASS:
            FILENGINE(pAi64, pSi64, pRi64, ai64, si64, int64);
            return;
        case mxUINT64_CLASS:
            FILENGINE(pAu64, pSu64, pRu64, au64, su64, uint64);
            return;
        case mxINT32_CLASS:
            FILENGINE(pAi32, pSi32, pRi32, ai32, si32, int32);
            return;
        case mxUINT32_CLASS:
            FILENGINE(pAu32, pSu32, pRu32, au32, su32, uint32);
            return;
        case mxCHAR_CLASS:
        case mxINT16_CLASS:
            FILENGINE(pAi16, pSi16, pRi16, ai16, si16, int16);
            return;
        case mxUINT16_CLASS:
            FILENGINE(pAu16, pSu16, pRu16, au16, su16, uint16);
            return;
        case mxINT8_CLASS:
            FILENGINE(pAi08, pSi08, pRi08, ai08, si08, int08);
            return;
        case mxLOGICAL_CLASS:
        case mxUINT8_CLASS:
            FILENGINE(pAu08, pSu08, pRu08, au08, su08, uint08);
            return;
    }
}