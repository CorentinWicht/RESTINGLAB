/**************************************************************************
 * MATLAB mex function mcolonmex
 *
 * [id1 id2 ...] = getclassid(A1, A2, ...)
 *
 * INPUTS: A1, A2, numerical arrays
 * OUTPUTS: id1, id2 ID of the class, started from 1, ordered in following:
 *          {'logical' 'char' ...
 *           'uint8' 'int8' 'uint16' 'int16' ...
 *           'uint32' 'int32' 'uint64' 'int64' ...
 *           'single' 'double'};
 *
 * Compile:
 *      mex getclassid.c on 32-bit
 *      mex -largeArrayDims getclassid.c on 64-bit
 *
 * MEX engine used by CASTARRAYS.M
 *
 * Author: Bruno Luong <brunoluong@yahoo.com>
 * Date: 31-Dec-2010
 *
**************************************************************************/

#include "mex.h"
        
#define R plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    mwSize i;
    double *pR;
            
    R = mxCreateNumericMatrix(1, nrhs, mxDOUBLE_CLASS, mxREAL);
    pR = mxGetPr(R);
    for (i=0; i<nrhs; i++) {
        switch (mxGetClassID(prhs[i])) {
            case mxLOGICAL_CLASS:
                pR[i] = 1;
                continue;
            case mxCHAR_CLASS:
                pR[i] = 2;
                continue;
            case mxUINT8_CLASS:
                pR[i] = 3;
                continue;
            case mxINT8_CLASS:
                pR[i] = 4;
                continue;
            case mxUINT16_CLASS:
                pR[i] = 5;
                continue;
            case mxINT16_CLASS:
                pR[i] = 6;
                continue;
            case mxUINT32_CLASS:
                pR[i] = 7;
                continue;
            case mxINT32_CLASS:
                pR[i] = 8;
                continue;
            case mxUINT64_CLASS:
                pR[i] = 9;
                continue;
            case mxINT64_CLASS:
                pR[i] = 10;
                continue;
            case mxSINGLE_CLASS:
                pR[i] = 11;
                continue;
            case mxDOUBLE_CLASS:
                pR[i] = 12;
                continue;
        }
    }
  
    return;
    
}
