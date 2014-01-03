/*
 VBviterbi_log.c
 =========================================================================
 
 Copyright (C) 2012 Martin Linden, E-mail: bmelinden@gmail.com

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

/*
 Matlab function definition:
 s=VBviterbi_log(lnQ,lnqst) 

% Run the Viterbi algorithm using the log of the transition matrix (lnQ)
% and emission likelihood lnqst, e.g., find the most likely path for the
% following distribution 
% p(s(1:T)) ~  exp( sum_t  lnqst(t,s(t)) + lnQ(s(t-1),s(t)) ).
%
% This is intended for use with variation treatments of Hidden Markov
% Models, where the log of the path distribution comes out naturally from
% the compuation. The point of using th elog representation is to avoid
% numerical under- or overflow for models where the entries in lnQ and
% lnqst have very large magnitude.
%
% Martin Lindén, 2012-01-12


 *
 */

#include <math.h>
#include "mex.h"
#include "matrix.h"

/* input Arguments */
#define	LNQ_IN prhs[0]
#define	LNQST_IN prhs[1]


#define	S_OUT plhs[0]

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]){
    double *lnQ,*lnqst; /* inputparameters */
    double *s;          /* output parameter */
    int T,N,t,j,k,ss;                       /* temporary variables */
    double Z;
    double *lnp0, *lnp1,*lnpp;
    double *maxPrev;
    mxArray *mxmaxPrev,*mxpp, *mxp1, *mxp0; /* temporary arrays */
    
    /* check number of input/output arguments */
    if (nlhs != 1)
        mexErrMsgTxt("One output argument required.");
    if (nrhs != 2)
        mexErrMsgTxt("Two input argument required.");
    
    /* size of input variables */
    T = mxGetM(LNQST_IN); /* number of rows = number of time points*/
    N = mxGetN(LNQST_IN); /* number of columns = number of states*/
    
    /* check that the other variables have consistent sizes */
    if ( mxGetM(LNQ_IN) != N || mxGetN(LNQ_IN) != N)
        mexErrMsgTxt("Q is not an N by N matrix (N = # columns in lnqst).");
    
    /* retrieve input data */
    lnQ=mxGetPr(LNQ_IN);
    lnqst=mxGetPr(LNQST_IN);
    
    /* Create an mxArray for the output data */
    S_OUT = mxCreateNumericMatrix(T,1,mxDOUBLE_CLASS, mxREAL);
    s=mxGetPr(S_OUT); /* pointer to output array*/
    
    /* allocate temporary arrays */
    mxp0=mxCreateDoubleMatrix(1,N, mxREAL);
    mxp1=mxCreateDoubleMatrix(1,N, mxREAL);
    lnp0=mxGetPr(mxp0);
    lnp1=mxGetPr(mxp1);
    
    mxmaxPrev=mxCreateNumericMatrix(T,N,mxDOUBLE_CLASS,mxREAL);
    maxPrev=mxGetPr(mxmaxPrev);
    mxpp=mxCreateDoubleMatrix(N,1, mxREAL);
    lnpp=mxGetPr(mxpp);
    
    /* start of actual algorithm */
    /*lnP1=lnqst(1,:)-mean(lnqst(1,:)); % initial probability, not normalized */
    Z=0.0;
    for(k=0;k<N;k++){
        Z=Z+lnqst[k*T]/N;
    }
    for(k=0;k<N;k++){
        lnp1[k]=lnqst[k*T]-Z;
    }
    for(t=1;t<T;t++){
        for(k=0;k<N;k++){/*lnP0=lnP1;lnP1=zeros(1,N);*/
            lnp0[k]=lnp1[k];
            lnp1[k]=0.0;
        }
        for(j=0;j<N;j++){
            for(k=0;k<N;k++){
                /* lnPP(kV)=lnP0(kV)+lnQ(kV,jV)+lnqst(tV,jV); 
                 *  % probability of most likely path that ends with kV -> jV*/
                lnpp[k]=lnp0[k]+lnQ[k+j*N]+lnqst[t+j*T];
            }
        	/* probability of previous state before ending up at jV.
            [lnP1(jV),         MaxPrev(tV,jV)]=max(lnPP); */
            lnp1[j]=lnpp[0];
            maxPrev[t+j*T]=0;
            for(k=1;k<N;k++){
                if(lnpp[k]>lnp1[j]){
                    lnp1[j]=lnpp[k];
                    maxPrev[t+j*T]=k;
                }
            }
        }
        
        Z=0.0;         /*lnP1=lnP1-mean(lnP1); % rescale*/
        for(k=0;k<N;k++){
            Z=Z+lnp1[k]/N;
        }        
        for(k=0;k<N;k++){
            lnp1[k]=lnp1[k]-Z;
        }
    }
    s[T-1]=0;     /* [~,S(T)]=max(lnP1);  */
    Z=lnp1[0];
    for(k=1;k<N;k++){
        if(lnp1[k]>Z){
            Z=lnp1[k];
            s[T-1]=k;
        }
    }
    /* for tV=T-1:-1:1
     * S(tV)=MaxPrev(tV+1,S(tV+1));
     * end*/
    for(t=T-2;t>=0;t--){
        ss=(int)s[t+1];
        s[t]=maxPrev[t+1+T*ss];
    }    
    /* finished actual algorithm */
    /* transform back to matlab's index convention */
    for(t=0;t<T;t++){
        s[t]=s[t]+1;
    }
    /* destroy temporary arrays */
    mxDestroyArray(mxp0);
    mxDestroyArray(mxp1);
    mxDestroyArray(mxpp);
    mxDestroyArray(mxmaxPrev);
}
