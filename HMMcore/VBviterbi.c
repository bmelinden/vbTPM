/*function s=VBviterbi(Q,qst) 
% s=VBviterbi(Q,qst)

 most likely trajectory by the Viterbi algorithm for a trajectory
  where Q is the transition matrix, ans qst is the emission
  likelihood. Constructed for use with ensemble learning of Hidden
  Markov models, M.L. 2011-12-20 

  This is the c-version of the m-fil sViterbi.m 
  compilation:
  mex('VBviterbi.c',['-' computer('arch')],'-O')
 
   M.L. 2011-12-21 cpoied, and added some robustness against pp(:)=0 some some t.

copyright notice
=========================================================================
VBviterbi.c
=========================================================================

Copyright (C) 2013 Martin Lindén

E-mail: bmelinden@gmail.com
=========================================================================
This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or any later
version.   
This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

Additional permission under GNU GPL version 3 section 7

If you modify this Program, or any covered work, by linking or combining it
with Matlab or any Matlab toolbox, the licensors of this Program grant you 
additional permission to convey the resulting work.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

 */

#include <math.h>
#include "mex.h"
#include "matrix.h"

/* input Arguments */
#define	Q_IN prhs[0]
#define	QST_IN prhs[1]

/* output arguments */
#define	S_OUT plhs[0]
#define	PT_OUT plhs[1]

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
  double *Q,*qst; /* input parameters */
  double *s;      /* output parameter */
  int T,N,t,j,k,ss;                       /* temporary variables */
  double Z,*pt,*pp;          
  double *maxPrevious;
  mxArray *mxmaxPrevious,*mxpp, *mxpt; /* temporary arrays */
  
   /* check number of input/output arguments */
 if (nlhs != 1) 
 mexErrMsgTxt("One output argument required.");  
 if (nrhs != 2)
 mexErrMsgTxt("Two input argument required.");  

 /* size of input variables */
 T = mxGetM(QST_IN); /* number of rows = number of time points*/
 N = mxGetN(QST_IN); /* number of columns = number of states*/
 
 /* check that the other variables have consistent sizes */
 if ( mxGetM(Q_IN) != N || mxGetN(Q_IN) != N)
   mexErrMsgTxt("Q is not an N by N matrix (N = # columns in qst).");
 
 /* retrieve input data */
 Q=mxGetPr(Q_IN);
 qst=mxGetPr(QST_IN);
 
 /* Create an mxArray for the output data */
 S_OUT = mxCreateNumericMatrix(T,1,mxDOUBLE_CLASS, mxREAL);
 s=mxGetPr(S_OUT); /* pointer to output array*/
 /*PT_OUT = mxCreateNumericMatrix(T,N,mxDOUBLE_CLASS, mxREAL);
   pt=mxGetPr(PT_OUT);*/

/* allocate temporary arrays */
 mxpt=mxCreateDoubleMatrix(T,N, mxREAL);
 pt=mxGetPr(mxpt);
 mxmaxPrevious=mxCreateNumericMatrix(T,N,mxDOUBLE_CLASS,mxREAL);
 maxPrevious=mxGetPr(mxmaxPrevious);
 mxpp=mxCreateDoubleMatrix(N,1, mxREAL);
 pp=mxGetPr(mxpp);

/* start of actual algorithm */
/* pt(1,:)=qst(1,:)/sum(qst(1,:)); % initial probability 
   pp=zeros(1,N); */
Z=0.0;
for(j=0;j<N;j++){
    Z=Z+qst[0+j*T];  
    pp[j*T]=0.0;
}
if(Z>0){
    for(j=0;j<N;j++){
        pt[j*T]=qst[0+j*T]/Z;
    }}else{
    mexEvalString("warning('VBviteri.c error: qst(1,:) not normalizable');");
    for(j=0;j<N;j++){
        pt[j*T]=1.0/N;
    }
    }

for(t=1;t<T;t++){
    for(j=0;j<N;j++){       /* current target state */
        pt[t+j*T]=0.0;
        maxPrevious[t+j*T]=maxPrevious[t-1+j*T]; /* default guess in case all pp=0 : same as previous */
        for(k=0;k<N;k++){   /* previous state */
            /* probability of most likely path that ends with k -> j : */
            pp[k]=pt[t-1+k*T]*Q[k+j*N]*qst[t+j*T]; 
            if (pp[k]>pt[t+j*T]){
                pt[t+j*T]=pp[k];
                maxPrevious[t+j*T]=k;
            }
        }
    }
    /* normalize pt[t,:] : */
    Z=0.0;
    for(j=0;j<N;j++){
        Z=Z+pt[t+j*T];
    }
    if(Z==0){
        mexPrintf("VBviterbi warning: pt(%d,:) not normalizable!\n",t);
        Z=1.0;
    }   
    for(j=0;j<N;j++){
        pt[t+j*T]=pt[t+j*T]/Z;
    }
 }
/* [~,s(T)]=max(pt(T,:));  */
Z=0.0;
for(k=0;k<N;k++){
    if(Z<pt[T-1+k*T]){
        Z=pt[T-1+k*T];
        s[T-1]=k;
}}
for(t=T-2;t>=0;t--){
    ss=(int)s[t+1];
    s[t]=maxPrevious[t+1+ss*T];
}
/* finished actual algorithm */
/* transform back to matlab's index convention */
for(t=0;t<T;t++){
  s[t]=s[t]+1;
 }
 
/* destroy temporary arrays */
mxDestroyArray(mxpt);
mxDestroyArray(mxpp);
mxDestroyArray(mxmaxPrevious);
}



















