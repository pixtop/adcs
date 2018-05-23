#include "mex.h"

void  subspace_iter1_(int*, int*, double*, double*, int*,
                            double*,
                            double*, double*,
                            int*,
                            double*, int*,
                            int* );


int ISEED[4] = {0,0,0,1};   /* initial seed for dlarnv() */
int IONE=1;

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  double *a, *V, *w, *res_ev, *z, *zi;
  int    *it_ev;
  double percentage, eps;
  int    ver, m, n, i, j, maxit, lwork, it, n_ev, ierr, mn;

  if(nrhs!=6) 
    mexErrMsgIdAndTxt( "MATLAB:subspace_iter_ev:invalidNumInputs",
                       "six inputs required.");
  if(nlhs!=4) 
    mexErrMsgIdAndTxt( "MATLAB:subspace_iter_ev:invalidNumOutputs",
                       "four outputs required.");

  /* get a pointer to the matrix */
  a = mxGetPr(prhs[0]);
  n = mxGetM(prhs[0]);
  
  /*  get the scalar input l */
  m = mxGetScalar(prhs[1]);
  mn = m*n;

  
  /*  get the scalar input p */
  ver = mxGetScalar(prhs[2]);

  /*  get the scalar input m */
  percentage = mxGetScalar(prhs[3]);

  /*  get the scalar input eps */
  eps = mxGetScalar(prhs[4]);

  /*  get the scalar input maxit */
  maxit = mxGetScalar(prhs[5]);
  
  V      = (double*)malloc(m*n*sizeof(double));
  w      = (double*)malloc(m*sizeof(double));
  res_ev = (double*)malloc(m*sizeof(double));
  it_ev  = (int*)malloc(m*sizeof(int));
  ISEED[0]=0; ISEED[1]=0; ISEED[2]=0; ISEED[3]=1; 
  dlarnv_(&IONE, ISEED, &mn, V);

  subspace_iter1_(&n, &m, a, &percentage, &maxit, &eps,
                  V, w, &n_ev, res_ev, it_ev,
                  &ierr );

  mexPrintf("---> %d %d  %d %d\n",n,m,n_ev, sizeof(int));

  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix( (mwSize)n, (mwSize)n_ev, mxREAL);
  /*  create a C pointer to a copy of the output matrix */
  z = mxGetPr(plhs[0]);

  for(j=0; j<n_ev; j++){
    for(i=0; i<n; i++){
      z[j*n+i] = V[j*n+i];
    }
  }


  /*  set the output pointer to the output matrix */
  plhs[1] = mxCreateDoubleMatrix( (mwSize)n_ev, (mwSize)1, mxREAL);
  /*  create a C pointer to a copy of the output matrix */
  z = mxGetPr(plhs[1]);

  for(j=0; j<n_ev; j++){
      z[j] = w[j];
  }

  /*  set the output pointer to the output matrix */
  plhs[2] = mxCreateDoubleMatrix( (mwSize)n_ev, (mwSize)1, mxREAL);
  /*  create a C pointer to a copy of the output matrix */
  z = mxGetPr(plhs[2]);

  for(j=0; j<n_ev; j++){
      z[j] = res_ev[j];
  }

    /*  set the output pointer to the output matrix */
  plhs[3] = mxCreateDoubleMatrix( (mwSize)n_ev, (mwSize)1, mxREAL);
  /*plhs[3] = mxCreateNumericMatrix((mwSize)n_ev, (mwSize) 1, mxINT32_CLASS, 0);*/
  /*  create a C pointer to a copy of the output matrix */
  zi = mxGetPr(plhs[3]);

  for(j=0; j<n_ev; j++){
      zi[j] = (double) it_ev[j];
  }

  free(V);
  free(w);
  free(res_ev);
  free(it_ev);

}
