**********************************************************************
|
| Program Name: EX1-GCMR-Sim1K-Sub100-Rho060.sas
|
| Program Version: 1
|
| Program Purpose: Fits Gauss Copula Marginal Regression to Normal
|                  and Binary data joinly. 
|
| SAS Version: 9.3
|
| Created By:     Thomas Drury
| Date:            August-2017
|
**********************************************************************;

ods html5 close;
ods listing gpath="Add Path to Store Plot Files";

libname brdata "Add Data Path" access = readonly;
libname results "Add Results Path";

%let num = 1;
%let dtype = rho60;

***********************************************************************;
*** FIT GAUSS COPULA MODEL                                          ***;
***********************************************************************;

proc sort data = brdata.ex1_sim1K_sub100_&dtype._&num.
          out  = indata&num.;
  by simulation subjid;
run;

ods select none;
ods results off;
ods graphics off;

proc mcmc data        = indata&num. 
          seed        = 12345678
          nbi         = 5000                        
          nmc         = 100000                      
          ntu         = 2000                                                     
          thin        = 10
          diagnostics = mcse
          outpost     = ex1_gcmr_sim1Ksub100&dtype._psam&num.
          statistics  = (summary interval)
          monitor     = (nt1 nt2 bt1 bt2 st1 st2 vt1 vt2
                         rt1 rt2 rho1 rho2 
                         td21 or21 dp21 p1 p2 br_pp1);

  by simulation;

  ods output   postsummaries = ex1_gcmr_sim1Ksub100&dtype._psum&num. 
               postintervals = ex1_gcmr_sim1Ksub100&dtype._pint&num.
               mcse          = ex1_gcmr_sim1Ksub100&dtype._mcse&num.;
 
  ****************************;
  *** MCMC PARAMETER SETUP ***;
  ****************************;
 
  array t[2,2] nt1 nt2
               bt1 bt2;
  array r[2] rt1 rt2;
  array s[2] st1 st2;

  parms bt1 0.1 bt2 0; 
  parms nt1 -50 nt2 -100;
  parms st1 50 st2 50;
  parms rt1 0 rt2 0;

  prior nt: ~ normal(0,sd=1000);
  prior bt: ~ normal(0,sd=1000);
  prior st: ~ igamma(shape=0.001,scale=0.001);
  prior rt: ~ uniform(-1,1);


  ****************************;
  ***        MODELS        ***;
  ****************************;

  *** MODEL FOR MEANS ***;
  xbeta1 = t[1,trtcd]; 
  xbeta2 = t[2,trtcd]; 

  *** MODELS FOR SD AND RHO ***;
  sd     = s[trtcd];
  theta  = r[trtcd];
  
  *** LINK FUNCTIONS ***;
  mu  = xbeta1;
  bp  = probnorm(xbeta2);                 


  ******************************;
  *** BUILD JOINT LIKELIHOOD ***;
  ******************************;

  fy1     = pdf("normal",y1,mu,sd); 
  cdf_y1  = cdf("normal",y1,mu,sd);    

  u1  = min(max(cdf_y1,1.0e-15),1-1.0e-15); *** BOUND CDF TO ENSURE BETWEEN 0 AND 1 ***; 
  u2  = min(max((1-bp),1.0e-15),1-1.0e-15); *** BOUND CDF TO ENSURE BETWEEN 0 AND 1 ***;

  q1  = probit(u1);   
  q2  = probit(u2);   

  c_u1u2 = probnorm( (q2 - (theta*q1)) / (sqrt(1-(theta**2))) );

  like = (y2=0)*fy1*c_u1u2 + (y2=1)*fy1*(1-c_u1u2);
  ll = log(like);

  model general(ll);


  *****************************************;
  ***      POSTERIOR PROCESSING         ***;
  *****************************************;

  beginnodata;

   *** BINARY PROBABILITIES ***;
   p1 = probnorm(bt1);
   p2 = probnorm(bt2);

   *** POINT BISERIAL ESTIMATES ***;
   rho1 = (rt1*pdf("normal",probit(p1),0,1))/sqrt(p1*(1-p1));
   rho2 = (rt2*pdf("normal",probit(p2),0,1))/sqrt(p2*(1-p2));

   *** NORMAL VARIANCES ***;
   vt1 = st1**2;
   vt2 = st2**2;

   *** TREATMENT COMPARISONS ***;
   td21 = nt2-nt1;  
   dp21 = p2-p1;
   or21 = (p2/(1-p2)) / (p1/(1-p1)); 
   
   *** POSTERIOR BR PROB ***;
   br_pp1 = (td21 > 100) * (dp21 < 0.3);
  
  
  endnodata;

run;
ods graphics off;
ods results on;
ods select all;


**********************************************************************;
***                    OUTPUT RESULTS DATASETS                     ***;
**********************************************************************;

data results.ex1_gcmr_sim1Ksub100&dtype._psam&num.;
  set ex1_gcmr_sim1Ksub100&dtype._psam&num.;
run;

data results.ex1_gcmr_sim1Ksub100&dtype._psum&num.;
  set ex1_gcmr_sim1Ksub100&dtype._psum&num.;
run;

data results.ex1_gcmr_sim1Ksub100&dtype._pint&num.;
  set ex1_gcmr_sim1Ksub100&dtype._pint&num.;
run;

data results.ex1_gcmr_sim1Ksub100&dtype._mcse&num.;
  set ex1_gcmr_sim1Ksub100&dtype._mcse&num.;
run;


**********************************************************************;
***                        CLEAN UP WORK AREA                      ***;
**********************************************************************;

proc datasets lib = work kill nolist;
quit;
run;
