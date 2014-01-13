/* [Za,alpha,Zb,beta,qt]=VB7_VBE_forwback(Q,qst)
% performs the forward and backward sweeps of the VBE step

copyright notice
=========================================================================
VBforwback.c, part of the vbTPM package
=========================================================================
Copyright (C) 2014 Martin Lindén, E-mail: bmelinden@gmail.com
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
=========================================================================

T=size(qst,1);
N=size(qst,2);

Za=zeros(T,1);
alpha=zeros(T,N);
alpha(1,:)=qst(1,:);
Za(1)=sum(alpha(1,:));
alpha(1,:)=alpha(1,:)/Za(1);
for t=2:T
    alpha(t,:)=(alpha(t-1,:)*Q).*qst(t,:);
    Za(t)=sum(alpha(t,:));
    alpha(t,:)=alpha(t,:)/Za(t);    
end % forward iterations

beta=zeros(T,N);
Zb=zeros(T,1); % only for debugging reasons...
Zb(T)=N;      % only for debugging reasons...
beta(T,:)=ones(1,N)/N;
QT=Q';
for t=T-1:-1:1
    beta(t,:)=(beta(t+1,:).*qst(t+1,:))*QT;
    Zb(t)=sum(beta(t,:));
    beta(t,:)=beta(t,:)/sum(beta(t,:));
end % backward iterations

 *
 *
 *
 *
 *
 *
 */
#include <math.h>
#include "mex.h"
#include "matrix.h"

/* Input Arguments */
#define	Q_IN prhs[0]
#define	QST_IN   prhs[1]

/* output arguments */
#define	ZA_OUT      plhs[0]
#define	ALPHA_OUT	plhs[1]
#define	ZB_OUT      plhs[2]
#define	BETA_OUT	plhs[3]
#define QT_OUT      plhs[4]

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[]){

  int N,T,j,k,t;
  double *qst,*Q;                  /* input parameters */
  double *Za,*alpha,*Zb,*beta,*qt; /* output parameters */
  mxArray *mxQt,*mx_qtZ;
  double  *Qt,*qtZ;
  

 /* check number of input/output arguments */
 if (nlhs != 5) 
 mexErrMsgTxt("Five output argument required.");  
 if (nrhs != 2)
 mexErrMsgTxt("Two input argument required.");  

 /* size of input variables */
 T = mxGetM(QST_IN); /* number of rows    */
 N = mxGetN(QST_IN); /* number of columns */

 /* check that the other variables have consistent sizes */
 if ( mxGetM(Q_IN) != N || mxGetN(Q_IN) != N)
   mexErrMsgTxt("Q is not an N by N matrix.");
   
 /* Create an mxArray for the output data */
 ZA_OUT = mxCreateDoubleMatrix(T, 1, mxREAL);
 Za=mxGetPr(ZA_OUT); /* create pointer to real data in output */
 ALPHA_OUT = mxCreateDoubleMatrix(T, N, mxREAL);
 alpha=mxGetPr(ALPHA_OUT); /* create pointer to real data in output */
 ZB_OUT = mxCreateDoubleMatrix(T, 1, mxREAL);
 Zb=mxGetPr(ZB_OUT); /* create pointer to real data in output */
 BETA_OUT = mxCreateDoubleMatrix(T, N, mxREAL);
 beta=mxGetPr(BETA_OUT); /* create pointer to real data in output */
 QT_OUT   = mxCreateDoubleMatrix(T, N, mxREAL);
 qt=mxGetPr(QT_OUT);

 /* retrieve input data */
 qst=mxGetPr(QST_IN);
 Q=mxGetPr(Q_IN);
 /* create temporary matrix: Q transpose */
 mxQt=mxCreateDoubleMatrix(N, N, mxREAL);
 Qt=mxGetPr(mxQt);
 mx_qtZ=mxCreateDoubleMatrix(T, 1, mxREAL);
 qtZ=mxGetPr(mx_qtZ);
 
 for(j=0;j<N;j++){
   for(k=0;k<N;k++){
     Qt[j+k*N]=Q[k+j*N];
   }
 }

for(j=0;j<N;j++){
  Za[0]=Za[0]+qst[j*T];
 }
for(j=0;j<N;j++){
  alpha[j*T]=qst[j*T]/Za[0];
 }
 /* real computation */
 for(t=1;t<T;t++){
   /* alpha(t,:)=(alpha(t-1,:)*Q).*qst(t,:);
      Za(t)=sum(alpha(t,:)); */
   for(j=0;j<N;j++){
     for(k=0;k<N;k++){
       alpha[t+j*T]=alpha[t+j*T]+alpha[t-1+k*T]*Q[k+j*N]*qst[t+j*T];
     }
     Za[t]=Za[t]+alpha[t+j*T];
   }
   /* alpha(t,:)=alpha(t,:)/Za(t); */
   for(j=0;j<N;j++){
     alpha[t+j*T]=alpha[t+j*T]/Za[t];
   }
 }

Zb[T-1]=N;
for(j=0;j<N;j++){
  beta[T-1+j*T]=1.0/N;
 }
 for(t=T-2;t>=0;t--){
   /* beta(t,:)=(beta(t+1,:).*qst(t+1,:))*QT;
      Zb(t)=sum(beta(t,:)); */
   for(j=0;j<N;j++){
     for(k=0;k<N;k++){
       beta[t+j*T]=beta[t+j*T]+(beta[t+1+k*T]*qst[t+1+k*T])*Qt[k+j*N];
     }
     Zb[t]=Zb[t]+beta[t+j*T];
   }
   /*     beta(t,:)=beta(t,:)/sum(beta(t,:)); */
   for(j=0;j<N;j++){
     beta[t+j*T]=beta[t+j*T]/Zb[t];
   }
 }
mxDestroyArray(mxQt);
/* qt=alpha.*beta;  % qt(t,j) = <(z(t)> = P(z(t)=j)
 * qtZ=sum(qt,2);
 * isNanInf=(~isempty(find(qtZ==0,1))); % check that qt is normalizable
 * if(isNanInf)
 * error('VB7_VBEMiter:qt_not_normalizeable','non-overlapping alpha(t,:) and beta(t,:) generated in VBE step.')
 * end
 * for j=1:NN
 * qt(1:end,j)=qt(1:end,j)./qtZ(1:end,1); % add realmin to avoid qt=0?
 * end  */
for(j=0;j<N;j++){
    for(t=0;t<T;t++){        
        qt[t+j*T]=alpha[t+j*T]*beta[t+j*T];
        qtZ[t]=qtZ[t]+qt[t+j*T];
    }
}
for(t=0;t<T;t++){
    if(qtZ[t]<=0.0){
        mexErrMsgTxt("VB7_VBEMiter:qt_not_normalizeable in VB7_VBE_forwback_4T.c");
    }
}
for(j=0;j<N;j++){
    for(t=0;t<T;t++){        
        qt[j*T+t]=qt[j*T+t]/qtZ[t];
    }
}
mxDestroyArray(mx_qtZ);
}
