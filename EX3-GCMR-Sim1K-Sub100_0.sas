**********************************************************************
|
| Program Name: EX3-GCMR-Sim1K-Sub100.sas
|
| Program Version: 1
|
| Program Purpose: Analyses Data For Example 3.
|
| SAS Version: 9.4
|
|
| Modified By:     Thomas Drury
| Date:            August 2017
|
**********************************************************************;

*********************************************************************;
***                     SET UP LIBNAMES                          ***;
*********************************************************************;

ods html5 close;
ods listing gpath = "Add Path to store Plot files";

libname brdata  "Add Path to the DATA";
libname results "Add Path to where RESULTS are to be stored";

%let num = 0;

***********************************************************************;
*** FIT GAUSS COPULA MODEL                                          ***;
***********************************************************************;

proc sort data = brdata.ex3_sim1K_sub100_&num.
          out  = indata&num.;
  by simulation usubjid;
run;

ods select none;
ods results off;
ods graphics off;

proc mcmc data = indata&num. 
          seed       = &num.12345  
          nbi        = 1000 
          nmc        = 5000 
          ntu        = 1000 
          thin       = 5 
          diagnostics = mcse dic 
          statistics  = (summary interval)
          outpost     = ex3_gcmr_sim1Ksub100_psam&num.
          monitor     =(nt1 nt2 st1 st2 
                       bt1 bt2 
                       nbt1 nbt2 dt1 dt2 
                       t1_o12 t1_o13 t1_o23_1 t1_o23
                       t2_o12 t2_o13 t2_o23_1 t2_o23
                       t1_bp t2_bp
                       t1_nbr t2_nbr   
                       );
    by simulation;
    
    ods output postsummaries = ex3_gcmr_sim1Ksub100_psum&num.;
    ods output postintervals = ex3_gcmr_sim1Ksub100_pint&num.;
    ods output mcse          = ex3_gcmr_sim1Ksub100_mcse&num.;
    ods output dic           = ex3_gcmr_sim1Ksub100_dic&num.;
    
        
    ****************************;
    *** MCMC PARAMETER SETUP ***;
    ****************************;
    
    array t[3, 2] nt1 nt2 
                  bt1 bt2 
                  nbt1 nbt2;
                  
    array d[2] dt1 dt2;
    
    array s[2] st1 st2;
    
    array o[3, 2] t1_o12 t2_o12 
                  t1_o13 t2_o13 
                  t1_o23_1 t2_o23_1;
    
    parms nt1 10 nt2 10 bt1 0.3 bt2 0.4 nbt1 1 nbt2 1;
    parms dt1 2 dt2 2;
    parms st1 50 st2 50;
    parms t1_o12 0 t2_o12 0 t1_o13 0 t2_o13 0 t1_o23_1 0 t2_o23_1 0;
    
    prior nt:  ~ normal(0, sd=1000);
    prior bt:  ~ normal(0, sd=1000);
    prior nbt: ~ normal(0, sd=1000);
    
    prior st: ~ igamma(shape=0.001, scale=0.001);
    prior dt: ~ uniform(0.1, 20);
    
    prior t1_o: ~ uniform(-1, 1);
    prior t2_o: ~ uniform(-1, 1);
    
    ****************************;
    ***        MODELS        ***;
    ****************************;
    
    *** MODEL FOR MEANS ***;
    xbeta1 = t[1, trtcd];
    xbeta2 = t[2, trtcd];
    xbeta3 = t[3, trtcd];
    
    *** MODELS FOR SD AND RHO ***;
    sd     = s[trtcd];
    disp   = d[trtcd];
    om12   = o[1, trtcd];
    om13   = o[2, trtcd];
    om23_1 = o[3, trtcd];

    *** LINK FUNCTIONS ***;
    mu   = xbeta1;
    bp   = probnorm(xbeta2);
    k    = 1 / disp;
    rate = exp(xbeta3);
    nbp  = k / (k + rate);

    ******************************;
    *** BUILD JOINT LIKELIHOOD ***;
    ******************************;

    ***************************;
    *** CONTINOUS VARIABLES ***;
    ***************************;

    pdf_y1=pdf("normal", y1, mu, sd);
    cdf_y1=cdf("normal", y1, mu, sd);

    **************************;
    *** DISCRETE VARIABLES ***;
    **************************;

    *** CALCULATE THE CDF FUNCTIONS FOR EACH VAR ***;
    cdf_y31=cdf("negbin", y3, nbp, k);
    cdf_y32=cdf("negbin", (y3-1), nbp, k);

    *** BOUND TO BE INSIDE LIMITS FOR PROBIT ***;

    u1  = min(max(cdf_y1, 1.0e-15), 1-1.0e-15);
    u21 = min(max((1-bp), 1.0e-15), 1-1.0e-15);
    u31 = min(max(cdf_y31, 1.0e-15), 1-1.0e-15);
    u32 = min(max(cdf_y32, 1.0e-15), 1-1.0e-15);

    *** CALCULATE COPULA TERMS CONDITIONAL ON Y1  ***;
    q21_1 = (probit(u21) - om12*probit(u1)) / (sqrt(1-(om12**2)));
    q31_1 = (probit(u31) - om13*probit(u1)) / (sqrt(1-(om13**2)));
    q32_1 = (probit(u32) - om13*probit(u1)) / (sqrt(1-(om13**2)));


    *** PAIRWISE DIFFERENCE THE DISCRETE VARIABLES ***;

    *** COMBINE INTO JOINT LIKELIHOOD ***;
    if y2=0 and y3=0 then do;
     c11_1 = probbnrm(q21_1, q31_1, om23_1);
     likelihood = pdf_y1*c11_1;
    end;
    else if y2=1 and y3=0 then do;
      u31_1 = probnorm(q31_1);
      c11_1 = probbnrm(q21_1, q31_1, om23_1);
      likelihood = pdf_y1*(u31_1 - c11_1);
    end;
    else if y2=0 and y3 gt 0 then do;
      c11_1 = probbnrm(q21_1, q31_1, om23_1);
      c21_1 = probbnrm(q21_1, q32_1, om23_1);
      likelihood = pdf_y1*(c11_1 - c21_1);
    end;
    else if y2=1 and y3 gt 0 then do;
      u31_1 = probnorm(q31_1);
      u32_1 = probnorm(q32_1);
      c11_1 = probbnrm(q21_1, q31_1, om23_1);
      c21_1 = probbnrm(q21_1, q32_1, om23_1);
      likelihood = pdf_y1*(u31_1 - u32_1 - c11_1 + c21_1);
    end;
    
    ll = log(likelihood);
    model general(ll);

    *****************************************;
    ***      POSTERIOR PROCESSING         ***;
    *****************************************;

    beginnodata;

    *** BACK CALCULATE T1_O23 AND T2_O23 FROM CONDITIONAL CORR FORMULA  ***;
    t1_o23 = (t1_o23_1 * sqrt(1-(t1_o13**2)) * sqrt(1-(t1_o12**2))) + (t1_o13*t1_o12);
    t2_o23 = (t2_o23_1 * sqrt(1-(t2_o13**2)) * sqrt(1-(t2_o12**2))) + (t2_o13*t2_o12);

    *** BACK CALCULATE BINARY PROBABILITY ***;
    t1_bp = probnorm(bt1);
    t2_bp = probnorm(bt2);
    
    *** BACK CALCULATE NEGBIN RATES ***;
    t1_nbr = exp(nbt1);
    t2_nbr = exp(nbt2);
        
    endnodata;
    
run;

ods select all;
ods results on;

**********************************************************************;
***                    OUTPUT RESULTS DATASETS                     ***;
**********************************************************************;

data results.ex3_gcmr_sim1Ksub100_psam&num.;
  set ex3_gcmr_sim1Ksub100_psam&num.;
run;

data results.ex3_gcmr_sim1Ksub100_psum&num.;
  set ex3_gcmr_sim1Ksub100_psum&num.;
run;

data results.ex3_gcmr_sim1Ksub100_pint&num.;
  set ex3_gcmr_sim1Ksub100_pint&num.;
run;

data results.ex3_gcmr_sim1Ksub100_mcse&num.;
  set ex3_gcmr_sim1Ksub100_mcse&num.;
run;

data results.ex3_gcmr_sim1Ksub100_dic&num.;
  set ex3_gcmr_sim1Ksub100_dic&num.;
run;


**********************************************************************;
***                        CLEAN UP WORK AREA                      ***;
**********************************************************************;

proc datasets lib = work kill nolist;
quit;
run;
