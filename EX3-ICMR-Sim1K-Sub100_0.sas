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
ods listing gpath = "Add Path to Store Plot Files";

libname brdata  "Add Data Path";
libname results "Add Results Path";

%let num = 0;

***********************************************************************;
*** FIT GAUSS COPULA MODEL                                          ***;
***********************************************************************;

proc sort data = brdata.ex3_sim1K_sub100_&num.
          out  = indata&num.;
  by simulation usubjid;
run;

ods select none;
ods results on;
ods graphics on;

proc mcmc data = indata&num. 
          seed       = &num.12345  
          nbi        = 1000 
          nmc        = 5000 
          ntu        = 1000 
          thin       = 5 
          diagnostics = mcse dic 
          statistics  = (summary interval)
          outpost     = ex3_icmr_sim1Ksub100_psam&num.
          monitor     =(nt1 nt2 st1 st2 
                       bt1 bt2 
                       nbt1 nbt2 dt1 dt2 
                       t1_bp t2_bp
                       t1_nbr t2_nbr   
                       );
                                              
    by simulation;
    
    ods output postsummaries = ex3_icmr_sim1Ksub100_psum&num.;
    ods output postintervals = ex3_icmr_sim1Ksub100_pint&num.;
    ods output mcse          = ex3_icmr_sim1Ksub100_mcse&num.;
    ods output dic           = ex3_icmr_sim1Ksub100_dic&num.;
    
        
    ****************************;
    *** MCMC PARAMETER SETUP ***;
    ****************************;
    
    array t[3, 2] nt1 nt2 
                  bt1 bt2 
                  nbt1 nbt2;
                  
    array d[2] dt1 dt2;
    
    array s[2] st1 st2;
    
    parms nt1 10 nt2 10 bt1 0.3 bt2 0.4 nbt1 1 nbt2 1;
    parms dt1 2 dt2 2;
    parms st1 50 st2 50;
    
    prior nt:  ~ normal(0, sd=1000);
    prior bt:  ~ normal(0, sd=1000);
    prior nbt: ~ normal(0, sd=1000);
    
    prior st: ~ igamma(shape=0.001, scale=0.001);
    prior dt: ~ uniform(0.1, 20);
    
    
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

    *** LINK FUNCTIONS ***;
    mu   = xbeta1;
    bp   = probnorm(xbeta2);
    k    = 1 / disp;
    rate = exp(xbeta3);
    nbp  = k / (k + rate);

    ******************************;
    *** BUILD JOINT LIKELIHOOD ***;
    ******************************;

    pdf_y1   = pdf("normal",y1,mu,sd); 
    pdf_y2   = pdf("bern",y2,bp); 
    pdf_y3   = pdf("negbin",y3,nbp,k); 
 
    ll = log(pdf_y1*pdf_y2*pdf_y3);
    model general(ll);
    
    
    *****************************************;
    ***      POSTERIOR PROCESSING         ***;
    *****************************************;

    beginnodata;

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

data results.ex3_icmr_sim1Ksub100_psam&num.;
  set ex3_icmr_sim1Ksub100_psam&num.;
run;

data results.ex3_icmr_sim1Ksub100_psum&num.;
  set ex3_icmr_sim1Ksub100_psum&num.;
run;

data results.ex3_icmr_sim1Ksub100_pint&num.;
  set ex3_icmr_sim1Ksub100_pint&num.;
run;

data results.ex3_icmr_sim1Ksub100_mcse&num.;
  set ex3_icmr_sim1Ksub100_mcse&num.;
run;

data results.ex3_icmr_sim1Ksub100_dic&num.;
  set ex3_icmr_sim1Ksub100_dic&num.;
run;


**********************************************************************;
***                        CLEAN UP WORK AREA                      ***;
**********************************************************************;

proc datasets lib = work kill nolist;
quit;
run;
