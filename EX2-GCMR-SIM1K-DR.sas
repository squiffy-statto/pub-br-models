**********************************************************************
|
| Program Name: EX2-GCMR-DR.sas
|
| Program Purpose: Estimates joint distribution of normal and binary 
|                  outcomes using Gaussian copulas for each simulated dataset
|                   
| SAS Version: 9.3
|
| Authors: Maria Costa and Thomas Drury
|
|********************************************************************;

*********************************************************************;
***                     SET UP LIBNAMES                          ***;
*********************************************************************;

ods html close;
ods listing gpath="add GPATH name";

libname brdata "add DATA library" access = readonly;
libname results "add RESULTS library";


***********************************************************************;
*** FIT GAUSS COPULA MODEL                                          ***;
***********************************************************************;

proc sort data = brdata.ex2_sim1k_sub50_dr
          out  = indata;
  by simulation subjid;
run;

ods select none;
ods results off;
ods graphics off;

proc mcmc data        = indata
          seed        = 12345678
          nbi         = 5000                        
          nmc         = 100000                      
          ntu         = 2000                                                     
          thin        = 10
          dic
          diagnostics = mcse
          statistics  = (summary interval)
          outpost = ex2_gcmr_sim1K_dr_psam
          monitor = (e_max e_d50
                     e_0 
                     sigma1-sigma6
                     alpha beta
                     theta1-theta6
                     p_d1-p_d6
                     rho_d1-rho_d6
                     m_d1-m_d6
                     variance1-variance6);
  by simulation;

  ods output postsummaries = ex2_gcmr_sim1K_dr_psum 
             postintervals = ex2_gcmr_sim1K_dr_pint
             mcse          = ex2_gcmr_sim1K_dr_mcse;
 
  ****************************;
  *** MCMC PARAMETER SETUP ***;
  ****************************;

  parms e_max 20 e_d50 1 e_0 -5;
  parms sigma1-sigma6 100;
  parms alpha 0 beta 1;
  parms theta1-theta6 0;

  prior e_0    ~ uniform(-1000,-1);
  prior e_d50  ~ uniform(0.001,8);
  prior e_max  ~ normal(0,sd=1000);                
  prior alpha  ~ normal(0,sd=1000);
  prior beta   ~ normal(0,sd=1000);
  prior sigma1-sigma6 ~ igamma(shape = 0.001, scale= 0.001);
  prior theta1-theta6 ~ uniform(-1,1);


  ****************************;
  ***        MODELS        ***;
  ****************************;

  *** MODEL FOR MEANS ***;
  lp1 = e_0 + ( (e_max*dose) / (e_d50+dose) ); 
  lp2 = alpha + (beta*dose); 

  *** MODEL FOR SIGMA ***;
  if dose = 0 then sdev = sigma1;
  else if dose = 0.3 then sdev = sigma2;
  else if dose = 0.7 then sdev = sigma3;
  else if dose = 1   then sdev = sigma4;
  else if dose = 4   then sdev = sigma5;
  else if dose = 6   then sdev = sigma6;

  *** MODEL FOR THETA ***;
  if dose = 0 then theta = theta1;
  else if dose = 0.3 then theta = theta2;
  else if dose = 0.7 then theta = theta3;
  else if dose = 1   then theta = theta4;
  else if dose = 4   then theta = theta5;
  else if dose = 6   then theta = theta6;
 
  *** LINK FUNCTIONS ***;
  mu  = lp1;
  bp  = probnorm(lp2);                 


  ******************************;
  *** BUILD JOINT LIKELIHOOD ***;
  ******************************;

  pdf_y1  = pdf("normal",y1,mu,sdev); 
  cdf_y1  = cdf("normal",y1,mu,sdev);    

  u1  = min(max(cdf_y1,1.0e-15),1-1.0e-15); *** BOUND CDF TO ENSURE BETWEEN 0 AND 1 ***; 
  u2  = min(max((1-bp),1.0e-15),1-1.0e-15); *** BOUND CDF TO ENSURE BETWEEN 0 AND 1 ***;

  q1  = probit(u1);   
  q2  = probit(u2);   

  c_u1u2 = probnorm( (q2 - (theta*q1)) / (sqrt(1-(theta**2))) );

  like = (y2=0)*pdf_y1*c_u1u2 + (y2=1)*pdf_y1*(1-c_u1u2);
  ll = log(like);

  model general(ll);

  beginnodata;


   *** BINARY PROBABILITIES ***;
   p_d1 = probnorm(alpha + beta*0);
   p_d2 = probnorm(alpha + beta*0.3);
   p_d3 = probnorm(alpha + beta*0.7);
   p_d4 = probnorm(alpha + beta*1);
   p_d5 = probnorm(alpha + beta*4);
   p_d6 = probnorm(alpha + beta*6);

   *** POINT BISERIAL (RHO) ESTIMATES ***;
   rho_d1 = (theta1*pdf("normal",probit(p_d1),0,1))/sqrt(p_d1*(1-p_d1));
   rho_d2 = (theta2*pdf("normal",probit(p_d2),0,1))/sqrt(p_d2*(1-p_d2));
   rho_d3 = (theta3*pdf("normal",probit(p_d3),0,1))/sqrt(p_d3*(1-p_d3));
   rho_d4 = (theta4*pdf("normal",probit(p_d4),0,1))/sqrt(p_d4*(1-p_d4));
   rho_d5 = (theta5*pdf("normal",probit(p_d5),0,1))/sqrt(p_d5*(1-p_d5));
   rho_d6 = (theta6*pdf("normal",probit(p_d6),0,1))/sqrt(p_d6*(1-p_d6));

   *** MEAN RESPONSE ESTIMATES ***;
   m_d1 = e_0 + ((e_max*0)/(e_d50+0));
   m_d2 = e_0 + ((e_max*0.3)/(e_d50+0.3));
   m_d3 = e_0 + ((e_max*0.7)/(e_d50+0.7));
   m_d4 = e_0 + ((e_max*1)/(e_d50+1));
   m_d5 = e_0 + ((e_max*4)/(e_d50+4));
   m_d6 = e_0 + ((e_max*6)/(e_d50+6));

   *** NORMAL VARIANCES ***;
   variance1 = sigma1**2;
   variance2 = sigma2**2;
   variance3 = sigma3**2;
   variance4 = sigma4**2;
   variance5 = sigma5**2;
   variance6 = sigma6**2;

  endnodata;



run;
ods graphics on;
ods results on;
ods select all;


**********************************************************************;
***                    OUTPUT RESULTS DATASETS                     ***;
**********************************************************************;

data results.ex2_gcmr_sim1K_dr1_psam ;
  set ex2_gcmr_sim1K_dr_psam ;
run;

data results.ex2_gcmr_sim1K_dr1_psum;
  set ex2_gcmr_sim1K_dr_psum;
run;

data results.ex2_gcmr_sim1K_dr1_pint;
  set ex2_gcmr_sim1K_dr_pint;
run;

data results.ex2_gcmr_sim1K_dr1_mcse;
  set ex2_gcmr_sim1K_dr_mcse;
run;

**********************************************************************;
***                        CLEAN UP WORK AREA                      ***;
**********************************************************************;

proc datasets lib = work kill nolist;
quit;
run;

