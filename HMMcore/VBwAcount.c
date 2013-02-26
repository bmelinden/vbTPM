/* wA=VB7_wAcount1m(alpha,beta,qst,Q)
 
performs the simple transition count in the VBE step (for the case
when there are no spurious states)

 */
#include <math.h>
#include "mex.h"
#include "matrix.h"

/* Input Arguments */
#define	ALPHA_IN prhs[0]
#define	BETA_IN	 prhs[1]
#define	QST_IN   prhs[2]
#define	Q_IN prhs[3]

/* output arguments */
#define	WA_OUT	plhs[0]


void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{

  int N,T,j,k,N2,T2,t;
  double *alpha,*beta,*qst,*Q; /* input parameters */
  double *wA;                  /* output parameter */
  double Z,*P;                    /* temporary variables */
  mxArray *mxP;
  
 /* check number of input/output arguments */
 if (nlhs != 1) 
 mexErrMsgTxt("One output argument required.");  
 if (nrhs != 4)
 mexErrMsgTxt("Four input argument required.");  

 /* size of input variables */
 T = mxGetM(ALPHA_IN); /* number of rows */
 N = mxGetN(ALPHA_IN); /* number of columns */

 /* check that the other variables have consistent sizes */
 if ( mxGetM(BETA_IN) != T || mxGetN(BETA_IN) != N)
   mexErrMsgTxt("Beta has different size than alpha.");
 if ( mxGetM(QST_IN) != T || mxGetN(QST_IN) != N)
   mexErrMsgTxt("qst has different size than alpha.");
 if ( mxGetM(Q_IN) != N || mxGetN(Q_IN) != N)
   mexErrMsgTxt("Q is not an N by N matrix.");
 
  
 /* Create an mxArray for the output data */
 WA_OUT = mxCreateDoubleMatrix(N, N, mxREAL);
 wA=mxGetPr(WA_OUT); /* create pointer to real data in output */
 
 /* retrieve input data */
 alpha=mxGetPr(ALPHA_IN);
 beta=mxGetPr(BETA_IN);
 qst=mxGetPr(QST_IN);
 Q=mxGetPr(Q_IN);
 
 /* create temporary matrix*/
 mxP=mxCreateDoubleMatrix(N, N, mxREAL);
 P=mxGetPr(mxP);
 
 /* real computation */
 for(t=1;t<T;t++){
     Z=0.0;
     for(j=0;j<N;j++){
         for(k=0;k<N;k++){
             P[j+k*N]=alpha[t-1+T*j]*Q[j+N*k]*qst[t+T*k]*beta[t+T*k];
             Z=Z+P[j+k*N];
         }
     }
     for(j=0;j<N*N;j++){
         wA[j]=wA[j]+P[j]/Z;
     }
 }
 mxDestroyArray(mxP);
}
