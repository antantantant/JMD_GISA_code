#include "mex.h"
#include "matrix.h"
#include <math.h>

#define dmod(x,y) ((x)-(y)*floor((x)/(y)))

static double xn = -1.;
double getran()
{
        static double a = 16807.;
        static double rm = 2147483647.; /* == 2^31-1 */
        static double seed = 123457.; /* good guess - or - 1234567*/

        if(xn== -1.) xn = dmod((a*seed),rm);
        else xn = dmod((a*xn),rm);
                                return((xn/(rm+1.)));
}
int reset_ran()
{
        xn=(-1);
        return(0);
}

void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
                 int nrhs, const mxArray *prhs[]) /* Input variables */
{   
    /* outputs: W */
    #define W plhs[0]

    /* inputs: s,p,mu,sigma,A,b,w0 */
    #define s prhs[0]
    #define p prhs[1]
    #define A prhs[2]
    #define b prhs[3]
    #define w0 prhs[4]
            
    int N, i, j, k, l, numConstraints, S, P;
    double *s_p, *p_p;
    double *W_p, *A_p, *w0_p, *boverA_p, *b_p, *as_p, *bounds_p;
    double lb, ub, lbb, ubb, random;
    mxArray *as, *boverA, *bounds, *input[1], *output[1];
    double a1, a2, a3, a4, a5, pp;
    double nlb, nub, t, y, prb;
    int sign,lbk, ubk;
    double A1, A2, A3, A4, A5, A6, B1, B2, B3, B4, B5, C1, C2, C3, C4, C5, C6, D1, D2, D3, D4;
    double p_low, p_high, q;

    
    
    S = (int)mxGetScalar(s);
    P = (int)mxGetScalar(p);
    
    N = 2*S;
    W = mxCreateDoubleMatrix(P, N, mxREAL);//fills up one column first
    W_p = mxGetPr(W);
    
    w0_p = mxGetPr(w0);
    for(i=0;i<P;i++){
        W_p[i] = w0_p[i];
    }
    

    A_p = (double *)mxGetPr(A);
    b_p = (double *)mxGetPr(b);
    
    numConstraints = (mwSize)mxGetN(A);
    as = mxCreateDoubleMatrix(numConstraints, 1, mxREAL);
    as_p = (double *)mxGetPr(as);
    
    boverA = mxCreateDoubleMatrix(numConstraints, 1, mxREAL);
    boverA_p = (double *)mxGetPr(boverA);
    bounds = mxCreateDoubleMatrix(2, 1, mxREAL);
    bounds_p = (double *)mxGetPr(bounds);
    
    // parameters for normal cdf
    
    a1 = 0.254829592;
    a2 = -0.284496736;
    a3 =  1.421413741;
    a4 = -1.453152027;
    a5 =  1.061405429;
    pp  =  0.3275911;

    
    // parameters for inverse normal cdf
    A1 = -3.969683028665376e+01;
    A2 = 2.209460984245205e+02;
    A3 = -2.759285104469687e+02;
    A4 = 1.383577518672690e+02;
    A5 = -3.066479806614716e+01;
    A6 = 2.506628277459239e+00;
    B1 = -5.447609879822406e+01;
    B2 = 1.615858368580409e+02;
    B3 = -1.556989798598866e+02;
    B4 = 6.680131188771972e+01;
    B5 = -1.328068155288572e+01;
    C1 = -7.784894002430293e-03;
    C2 = -3.223964580411365e-01;
    C3 = -2.400758277161838e+00;
    C4 = -2.549732539343734e+00;
    C5 = 4.374664141464968e+00;
    C6 = 2.938163982698783e+00;
    D1 = 7.784695709041462e-03;
    D2 = 3.224671290700398e-01;
    D3 = 2.445134137142996e+00;
    D4 = 3.754408661907416e+00;
    p_low = 0.02425;
    p_high = 1.-p_low;
    
    reset_ran();
    for(i = 1; i<N; i++){
//         if(dmod(i,1000)==0){mexPrintf("*");}
        for(j=0; j<P; j++){
            lb = -1000.;
            ub = 1000.;
            for(k=0; k<numConstraints; k++){
                as_p[k] = 0.0;
                for(l=0; l<j; l++){
                    as_p[k] += A_p[k*P+l] * W_p[i*P+l];
//                     mexPrintf("*j=%d,k=%d,l=%d,A=%g,W=%g,as=%g*",j,k,l,A_p[k*P+l],W_p[i*P+l],as_p[k]);
                }
                for(l=j+1;l<P;l++){
                    as_p[k] += A_p[k*P+l] * W_p[(i-1)*P+l];
//                     mexPrintf("*j=%d,k=%d,l=%d,A=%g,W=%g,as=%g*",j,k,l,A_p[k*P+l],W_p[(i-1)*P+l],as_p[k]);
                }
                boverA_p[k] = (b_p[k] - as_p[k])/A_p[k*P+j];
                
                if (A_p[k*P+j]>1e-3 && boverA_p[k]>lb){
                    lb = boverA_p[k];
//                     mexPrintf("|lb=%g,k=%d|",lb,k);
                    lbk = k;
                }
                else if (A_p[k*P+j]<-1e-3 && boverA_p[k]<ub){
                    ub = boverA_p[k];
//                     mexPrintf("|ub=%g,k=%d|",ub,k);
                    ubk = k;
                }
            }
//             if(i==1&&j==20){mexPrintf("|lb,lbk,ub,ubk,%g,%d,%g,%d|",lb,lbk,ub,ubk);}
            bounds_p[0] = lb;
            bounds_p[1] = ub;
            input[0] = bounds;
//             mexCallMATLAB(1, output, 1, input, "tmvn");
//             random = (double)mxGetScalar(output[0]);

            if(lb>ub){
//                mexPrintf("|p=%d,ub=%g,lb=%g|",j,ub,lb);
                random = (ub+lb)/2.;
            }
            else{
            // try truncated multivariate normal
                // Save the sign of x
                sign = 1.;
                if (lb < 0)
                    sign = -1.;
//                 mexPrintf("|lb=%g|",fabs(lb));
                lbb = fabs(lb)/sqrt(2.0);
                // A&S formula 7.1.26
                t = 1.0/(1.0 + pp*lbb);
                y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-lbb*lbb);
                
                nlb = 0.5*(1.0 + sign*y);
//                 mexPrintf("|nlb=%g|",nlb);

                sign = 1.;
                if (ub < 0)
                    sign = -1.;
                ubb = fabs(ub)/sqrt(2.0);
                // A&S formula 7.1.26
                t = 1.0/(1.0 + pp*ubb);
                y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-ubb*ubb);
//                 mexPrintf("|ub=%g, t=%g, y=%g|",ub,t,y);
                nub = 0.5*(1.0 + sign*y);
//                 mexPrintf("|nub=%g|",nub);
                
                prb = getran();
//                 mexPrintf("|prb=%g|",prb);
                prb = prb*(nub-nlb)+nlb;
//                 mexPrintf("|prb=%g|",prb);
                // inverse normal
                
                if (prb<1e-6){ // in case both ub and lb are on the tail
                    random = (ub+lb)/2.0;
//                     if(i==1&&j==20){mexPrintf("|prb,ub,lb=%g,ub,lb|",prb,ub,lb);}
                }
                else if (prb>(1-1e-6)){ // in case both ub and lb are on the tail
                    random = (ub+lb)/2.0;
//                     if(i==1&&j==20){mexPrintf("|prb,ub,lb=%g,ub,lb|",prb,ub,lb);}
                }    
                else if (prb<p_low){
                    q = sqrt(-2*log(prb));
                    random = (((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6)/((((D1*q+D2)*q+D3)*q+D4)*q+1);
//                     mexPrintf("|low:random=%g|",random);
                }
                else if (prb<p_high){
                    q = (prb - 0.5)*(prb - 0.5);
                    random = (((((A1*q+A2)*q+A3)*q+A4)*q+A5)*q+A6)*(prb-0.5)/(((((B1*q+B2)*q+B3)*q+B4)*q+B5)*q+1);
//                     mexPrintf("|med:random=%g|",random);
                }
                else if (prb<1){
                    q = sqrt(-2*log(1-prb));
                    random = -(((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6)/((((D1*q+D2)*q+D3)*q+D4)*q+1);
//                     mexPrintf("|high:random=%g|",random);
                }
            // try truncated multivariate normal
            }
//             if(i==1&&j==20){mexPrintf("|output:%g|",random);}

            W_p[i*P+j] = random;
//             W_p[i*P+j] = (lb+ub)/2;
//             mxDestroyArray(random);
        }
    }
    return;
}

