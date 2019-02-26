**********************************************************************
|
| Program Name: EX1-SRE-Sim1K-Sub100-Rho_1.sas
|
| Program Version: 1
|
| Program Purpose: Fits Shared Random Effect GLMM to model Normal
|                  and Binary data joinly. 
|
| SAS Version: 9.3
|
| Created By:      Maria Costa  
| Date:            April-2016
|
| Modified By:     Thomas Drury
| Date:            August-2017
|
**********************************************************************;

ods html5 close;
ods listing gpath="add GPATH line";

libname brdata "add DATA library" access = readonly;
libname results "add RESULTS library";

%let nsims = 1K;
%let dtype = rho20;  /* States which correlation scenario is being investigated */
%let npart = 1;

***********************************************************************;
*** FIT SHARED RANDOM EFFECTS GLMM                                  ***;
***********************************************************************;
options nosource nonotes;

proc sort data = brdata.ex1_sim&nsims._sub100_&dtype._&npart.
          out  = indata&npart.;
  by simulation subjid;
run;


ods results off;
ods select none;
ods graphics off;

proc mcmc data        = indata&npart. 
          nbi         = 5000                        
          nmc         = 250000                    
          ntu         = 3000                                                     
          thin        = 25
          seed        = 12345678
		  outpost     = ex1_sre_sim&nsims.sub100&dtype._psam&npart.
          diagnostics = mcse
          statistics  = (summary interval)
          monitor     = (n1 n2                /*** NORMAL MODEL ESTIMATES           ***/
                         b1m b2m              /*** BINARY MODEL ESTIMATES           ***/
                         t1var t2var          /*** VARIANCE OF NORMAL DATA FOR TRT1 AND TRT2 ***/
                         t1uvar t2uvar        /*** VARIANCE OF SHARED RE TRT1 AND TRT2 ***/      
                         rho1 rho2            /*** ESTIMATE OF CORR(Y1,Y2) FOR SUBJECTS ON TRT1 AND TRT2 ***/
                         totalvar1 totalvar2  /*** ESTIMATE OF TOTAL VARIANCE FOR SUBJECTS ON TRT1 AND TRT2 ***/
                         t1p t2p              /*** ESTIMATE OF P(Y=1|TRT=1) AND P(Y=1|TRT=2)  ***/
                         td or dp             /*** TREATMENT DIFFERENCE ODDS RATIO AND PROB DIFF ***/             
                         br_pp1);         

  by simulation;

  ods output   postsummaries = ex1_sre_sim&nsims.sub100&dtype._psum&npart. 
               postintervals = ex1_sre_sim&nsims.sub100&dtype._pint&npart.
               mcse          = ex1_sre_sim&nsims.sub100&dtype._mcse&npart.;

  array lp[2];
  array variances[2] t1var t2var;
  array u_variances[2] t1uvar t2uvar;

  parms n1 0 n2 0;
  parms b1 0 b2 0;
  parms t1var 100 t2var 100;
  parms t1uvar 1 t2uvar 1;
  
  prior n:  ~ normal(0,sd=1000);
  prior b:  ~ normal(0,sd=1000);
  prior t:  ~ igamma(shape = 0.001, scale = 0.001);
  
  lp[1] = n1*(trtcd=1) + n2*(trtcd=2) + u*sqrt(variances[trtcd]);   *** RESCALE U TO NORMAL SCALE ***;
  lp[2] = b1*(trtcd=1) + b2*(trtcd=2) + u;

  mu = lp[1];                      *** IDENTITY LINK FOR NORMAL MEAN ***;
  p  = probnorm(lp[2]);            *** REVERSE OF PROBIT LINK FOR BINARY PROBABILITY  ***;
 
  random u ~ normal(0,var=u_variances[trtcd]) subject = subjid;

  model y1 ~ normal(mu,var=variances[trtcd]);
  model y2 ~ binary(p);

  beginnodata;
   b1m    = b1/sqrt(1+t1uvar);    *** CORRECTION FACTOR TO CREATE MARGINAL ESTIMATE ***;
   b2m    = b2/sqrt(1+t2uvar);    *** CORRECTION FACTOR TO CREATE MARGINAL ESTIMATE ***;
   t1p    = probnorm(b1m);
   t2p    = probnorm(b2m);
   td     = n2-n1;
   dp     = t2p-t1p;
   or     = (t2p/(1-t2p)) / (t1p/(1-t1p)); 
   rho1   = (t1uvar/(1+t1uvar)) * ( pdf("normal",b1m,0,1) / (sqrt(probnorm(b1m)*(1-probnorm(b1m))) ));
   rho2   = (t2uvar/(1+t2uvar)) * ( pdf("normal",b2m,0,1) / (sqrt(probnorm(b2m)*(1-probnorm(b2m))) ));
   totalvar1 = t1var * t1uvar + t1var;
   totalvar2 = t2var * t2uvar + t2var;
   
   *** POSTERIOR BR PROB ***;
   br_pp1 = (td > 100) * (dp < 0.3);
   
  endnodata;

run;

ods graphics off;
ods results on;
ods select all;


**********************************************************************;
***                    OUTPUT RESULTS DATASETS                     ***;
**********************************************************************;
data results.ex1_sre_sim&nsims.sub100&dtype._psam&npart.;
  set ex1_sre_sim&nsims.sub100&dtype._psam&npart.;
run;

data results.ex1_sre_sim&nsims.sub100&dtype._psum&npart.;
  set ex1_sre_sim&nsims.sub100&dtype._psum&npart.;
run;

data results.ex1_sre_sim&nsims.sub100&dtype._pint&npart.;
  set ex1_sre_sim&nsims.sub100&dtype._pint&npart.;
run;

data results.ex1_sre_sim&nsims.sub100&dtype._mcse&npart.;
  set ex1_sre_sim&nsims.sub100&dtype._mcse&npart.;
run;

options source notes;

**********************************************************************;
***                        CLEAN UP WORK AREA                      ***;
**********************************************************************;

proc datasets lib = work kill nolist;
quit;
run;
