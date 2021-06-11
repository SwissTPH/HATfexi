#include <Rcpp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <utility>

using namespace Rcpp;
//using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// using STLs improves the Gillespie implementation. The "artisanal" version I had previously crashes R
// Here I count annual cummulative values of I1L and I1H (true incidence) as well as nb of treated people per year (at least
//in this version I don't track whether the treated cases come from AS or PD

//add omega_pd
  
// [[Rcpp::export]]
List hat_gillespie_fexi(List params){
  // remember that now tdisc=0 means 2013 (jut after data and detemrinistic dynamics finishes and start the stochastic dynamics)
  // So now for any tdisc, use the mean of past (regular, half or double) AS rate
  // for now, all rates are constant because are projections
  double tole=0.00001; // tolerance to count and save data up to the last year
  NumericVector r_ASvec= params["rateAS"];//a scalar repeated eahc year for projections beyond the data period
  //NumericVector r_ASvec = rAS;
  NumericVector r_pd1vec = params["r_pd1"]; //PD detection (and treatment) rate for stage 1, current first-line treatment
  NumericVector r_pd2vec = params["r_pd2"]; //PD detection (and treatment) rate for stage 2, current first-line treatment
  //for use of fexinidazole in passive detection:
  double psi = params["psi"]; // to account for the  lower compliance to fexi.
  // double r_pd1Fexi= omega_pd*r_pd1;
  // double r_pd2Fexi= omega_pd*r_pd2;
  // double r_asFexi= omega_as*rAS; // I assume r_as1Fexi=r_as2Fexi, I prefer define r_asFexi to avoid confusion
  // double e1= params["e1"]; // to consider in case great difference in the efficacy is observed from data; currently will explore 0.87
  // double e2= params["e2"];
  double phi1 = params["phi1"];
  double phi2 = params["phi2"];
  // Pull out current fixed params for easy reading
  double alpha = params["alpha"];
  double b = params["b"];
  double c_aL = params["c_aL"];
  double c_aH = params["c_aH"];
  double delta = params["delta"];
  double delta_a = params["delta_a"];
  // double epsilon = params["epsilon"];
  double eta = params["eta"];
  double f = params["f"];
  double gamma = params["gamma"];
  double gamma_aL = params["gamma_aL"];
  double gamma_aH = params["gamma_aH"];
  double mu = params["mu"];
  double mu_aL = params["mu_aL"];
  double mu_aH = params["mu_aH"];
  double mu_t = params["mu_t"];
  double mu_v = params["mu_v"];
  // double Nh0 = params["Nh0"];
  double nu = params["nu"];
  double sensitivity = params["sensitivity"];
  double sigma = params["sigma"];
  double sigma_aL = params["sigma_aL"];
  double sigma_aH = params["sigma_aH"];
  double specificity = params["specificity"];
  double xi = params["xi"];
  
  
  // Pull out current SAMPLED params for easy reading (i.e. Pars que sampleo al arrancar la deterministic part)
  // double AH1 = params["AH1"];
  // double AH2 = params["AH2"];
  double c_h = params["c_h"];
  // double coef2 = params["coef2"];
  double mu_gamma = params["mu_gamma"];
  // double N2N1 = params["N2N1"];
  // double VH1 = params["VH1"];
  // double VH2 = params["VH2"];
  
  int maxyears = as<int>(params["maxyears"]);  // SC: after maxyears, stop the simulation
  
  // chained operations are tricky in cpp
  // pull out list w/in list into its own object
  List initHuman = params["initHuman"]; 
  List initVector = params["initVector"];
  List initAnimal = params["initAnimal"];
  
  
// Define and initialise state variables and other variables that I will track (e.g. reported cases, treated cases)
  //Humans low-risk setting
  double S_L = initHuman["S_L"]; 
  double E_L = initHuman["E_L"];
  double I1_L = initHuman["I1_L"];
  double I2_L = initHuman["I2_L"];
  double T_L = initHuman["T_L"];
  int byYear_Incid_1L = initHuman["Incidence_1L"]; //initialise to inheritated integer values
  //int byYear_Active1 = initHuman["Active1"]; //idem
  //int byYear_Active2 = initHuman["Active2"]; //idem
  //int byYear_Passive_1L = initHuman["Passive_1L"]; //idem
  //int byYear_Passive_2L = initHuman["Passive_2L"]; //idem
  
  // int byYearReported_AS = byYear_Active1+byYear_Active1; // reported cases AS  ###### DEBERIAN SER NUMEROS ENTEROS!
  //int byYearReported_PD_LR = r_pd1*I1_L+r_pd2*I2_L;  //reported cases PD LR ###### DEBERIAN SER NUMEROS ENTEROS!
  
  //Humans high-risk setting
  double S_H = initHuman["S_H"];
  double E_H = initHuman["E_H"];
  double I1_H = initHuman["I1_H"];
  double I2_H = initHuman["I2_H"];
  double T_H = initHuman["T_H"];
  int byYear_Incid_1H = initHuman["Incidence_1H"];  //initialise to inheritated integer values
  //int byYear_Passive_1H = initHuman["Passive_1H"]; //idem
  //int byYear_Passive_2H = initHuman["Passive_2H"]; //idem
  //int byYearReported_PD_HR = r_pd1*I1_H+r_pd2*I2_H;  //reported cases PD LR ###### DEBERIAN SER NUMEROS ENTEROS!
  
  //Define "variables" to track reported cases by group, stage and year
  
  int cumI_1L = I1_L; 
  //int cumI_2L = I2_L;
  int cumI_1H = I1_H;
  //int cumI_2H = I2_H;
  
  // int cumReported_Active1 = byYear_Passive_1L;
  // int cumReported_Active2 = byYear_Passive_2L;  //this only applies for LR setting
  // int byYearReported_Passive1L = byYear_Passive_1L; //up to just before starting stochastic, these (reported and treated) match, so it makes sense
  // int byYearReported_Passive2L = byYear_Passive_2L; 
  // int byYearReported_Passive1H = byYear_Passive_1H; 
  // int byYearReported_Passive2H = byYear_Passive_2H; 
  
  // Define variables to store info on diagnosis
  int byYear_Diagnosis_1L = 0;
  int byYear_Diagnosis_2L = 0;
  int byYear_Diagnosis_1H = 0; 
  int byYear_Diagnosis_2H = 0;  
  
  
  // Define variables to store info on treatment outcomes
  int byYear_Active1_fexi = 0;
  int byYear_Active1_curr = initHuman["Active1"]; 
  int byYear_Active2_fexi = 0; 
  int byYear_Active2_curr = initHuman["Active2"]; 
  
  int byYear_Passive1_fexi = 0; 
  int byYear_Passive1_curr = initHuman["Passive_1"]; 
  int byYear_Passive2_fexi = 0;
  int byYear_Passive2_curr = initHuman["Passive_2"];
  
  //int byYear_Passive1L_fexi = 0; 
  //int byYear_Passive1L_curr = initHuman["Passive_1L"]; 
  //int byYear_Passive2L_fexi = 0;
  //int byYear_Passive2L_curr = initHuman["Passive_2L"];
  
  //int byYear_Passive1H_fexi = 0;
  //int byYear_Passive1H_curr = initHuman["Passive_1H"];
  //int byYear_Passive2H_fexi = 0;
  //int byYear_Passive2H_curr = initHuman["Passive_2H"];
  
  //Vectors
  double SvL = initVector["SvL"];
  double EvL = initVector["EvL"];
  double IvL = initVector["IvL"];
  double UvL = initVector["UvL"];
  
  double SvH = initVector["SvH"];
  double EvH = initVector["EvH"];
  double IvH = initVector["IvH"];
  double UvH = initVector["UvH"];
  
  //Animals
  double SaL = initAnimal["SaL"];
  double EaL = initAnimal["EaL"];
  double IaL = initAnimal["IaL"];
  double RaL = initAnimal["RaL"];
  
  double SaH = initAnimal["SaH"];
  double EaH = initAnimal["EaH"];
  double IaH = initAnimal["IaH"];
  double RaH = initAnimal["RaH"];
  
  // Initial values for vectors to store different aspects of detecting and treating with different drugs,
  // They store annual values (saved and reset to zero after each year completed) instead of cumulatives 
  // REMOVI TODO DE ESTA SUBSECCION
  
  // Initialise vectors to store every state variable
  std::vector<double> tt;
  
/*   // Human population
  std::vector<double> SS_L;
  std::vector<double> EE_L;
  std::vector<double> II1_L;
  std::vector<double> II2_L;
  std::vector<double> TT_L;
  
  std::vector<double> SS_H;
  std::vector<double> EE_H;
  std::vector<double> II1_H;
  std::vector<double> II2_H;
  std::vector<double> TT_H;
  
  // Vector population
  std::vector<double> SSvL;
  std::vector<double> EEvL;
  std::vector<double> IIvL;
  std::vector<double> UUvL;
  
  std::vector<double> SSvH;
  std::vector<double> EEvH;
  std::vector<double> IIvH;
  std::vector<double> UUvH;
  
  // Animal population
  std::vector<double> SSaL;
  std::vector<double> EEaL;
  std::vector<double> IIaL;
  std::vector<double> RRaL;
  
  std::vector<double> SSaH;
  std::vector<double> EEaH;
  std::vector<double> IIaH;
  std::vector<double> RRaH; */
  
// this to track values from treatment I am interested on
  std::vector<double> bbyYear_Incid_1L;
  std::vector<double> bbyYear_Incid_1H;
  
  std::vector<double> bbyYear_Diagnosis_1L;
  std::vector<double> bbyYear_Diagnosis_2L;  
  std::vector<double> bbyYear_Diagnosis_1H;
  std::vector<double> bbyYear_Diagnosis_2H;  
  
  std::vector<double> bbyYear_Active1_fexi;
  std::vector<double> bbyYear_Active1_curr;
  std::vector<double> bbyYear_Active2_fexi;
  std::vector<double> bbyYear_Active2_curr;
  
  std::vector<double> bbyYear_Passive1_fexi;
  std::vector<double> bbyYear_Passive1_curr;
  std::vector<double> bbyYear_Passive2_fexi;
  std::vector<double> bbyYear_Passive2_curr;  
  
  //std::vector<double> bbyYear_Passive1L_fexi;
  //std::vector<double> bbyYear_Passive1L_curr;
  //std::vector<double> bbyYear_Passive2L_fexi;
  //std::vector<double> bbyYear_Passive2L_curr;
  //std::vector<double> bbyYear_Passive1H_fexi;
  //std::vector<double> bbyYear_Passive1H_curr;
  //std::vector<double> bbyYear_Passive2H_fexi;
  //std::vector<double> bbyYear_Passive2H_curr;

  // std::vector<double> ccumReported_Active1; //vector to store annual values
  // std::vector<double> ccumReported_Active2;
  // std::vector<double> bbyYearReported_Passive1L;
  // std::vector<double> bbyYearReported_Passive2L;
  // std::vector<double> bbyYearReported_Passive1H;
  // std::vector<double> bbyYearReported_Passive2H;

      
      //Initialise time and other params
  double t = 2000;
  double tstart = 2000; //fixed start year to allow calculation of index to draw from rate vectors
  double tdisc = 2000; //indicate (approx) when store updated states
  //double tdisc2 = 2000; //r_AS at beginning of each month
  double tdisc3 = 2021; //to model the current treatments from 2000 to 2021
  //double onemonth=0.08333333333;
  
  do{
    // simulate vector control
    //if(t>2){mu_v=mu_v/0.99;}
	 
    if (t>=tdisc) {

      // add updated state values
      tt.push_back(t);

      /* //humans
      SS_L.push_back(S_L);
      EE_L.push_back(E_L);
      II1_L.push_back(I1_L);
      II2_L.push_back(I2_L);
      TT_L.push_back(T_L);

      SS_H.push_back(S_H);
      EE_H.push_back(E_H);
      II1_H.push_back(I1_H);
      II2_H.push_back(I2_H);
      TT_H.push_back(T_H);

      //Vectors
      SSvL.push_back(SvL);
      EEvL.push_back(EvL);
      IIvL.push_back(IvL);
      UUvL.push_back(UvL);

      SSvH.push_back(SvH);
      EEvH.push_back(EvH);
      IIvH.push_back(IvH);
      UUvH.push_back(UvH);

      //Animals
      SSaL.push_back(SaL);
      EEaL.push_back(EaL);
      IIaL.push_back(IaL);
      RRaL.push_back(RaL);

      SSaH.push_back(SaH);
      EEaH.push_back(EaH);
      IIaH.push_back(IaH);
      RRaH.push_back(RaH); */

      // Reported cases when using fexi (remember that also using current for some special groups)
      // store annual values
      bbyYear_Incid_1L.push_back(byYear_Incid_1L);
      bbyYear_Incid_1H.push_back(byYear_Incid_1H);
	  
      bbyYear_Diagnosis_1L.push_back(byYear_Diagnosis_1L);
	  bbyYear_Diagnosis_2L.push_back(byYear_Diagnosis_2L);
      bbyYear_Diagnosis_1H.push_back(byYear_Diagnosis_1H);
	  bbyYear_Diagnosis_2H.push_back(byYear_Diagnosis_2H);
      
      bbyYear_Active1_fexi.push_back(byYear_Active1_fexi);
      bbyYear_Active1_curr.push_back(byYear_Active1_curr);
      bbyYear_Active2_fexi.push_back(byYear_Active2_fexi);
      bbyYear_Active2_curr.push_back(byYear_Active2_curr);

      bbyYear_Passive1_fexi.push_back(byYear_Passive1_fexi);
      bbyYear_Passive1_curr.push_back(byYear_Passive1_curr);
      bbyYear_Passive2_fexi.push_back(byYear_Passive2_fexi);
      bbyYear_Passive2_curr.push_back(byYear_Passive2_curr);
      
      //bbyYear_Passive1L_fexi.push_back(byYear_Passive1L_fexi);
      //bbyYear_Passive1L_curr.push_back(byYear_Passive1L_curr);
      
      //bbyYear_Passive2L_fexi.push_back(byYear_Passive2L_fexi);
      //bbyYear_Passive2L_curr.push_back(byYear_Passive2L_curr);
      
      //bbyYear_Passive1H_fexi.push_back(byYear_Passive1H_fexi);
      //bbyYear_Passive1H_curr.push_back(byYear_Passive1H_curr);
      
      //bbyYear_Passive2H_fexi.push_back(byYear_Passive2H_fexi);
      //bbyYear_Passive2H_curr.push_back(byYear_Passive2H_curr);
      
      // byYearReported_Passive1L = r_pd1*I1_L;  //puede que no de entero !! el 'cum' es por anio !
      // byYearReported_Passive2L = r_pd2*cumI_2L;
      // byYearReported_Passive1H = r_pd1*cumI_1H;
      // byYearReported_Passive2H = r_pd2*cumI_2H;
      
      //std::cout <<  byYearReported_Passive1L << std::endl; // printing  
      //double aa = byYear_Passive1L_fexi + byYear_Passive1L_curr
      //std::cout <<  byYear_Passive1L_fexi << std::endl; // printing  
      //std::cout <<  byYear_Passive1L_curr << std::endl; // printing  
      
      //ccumReported_Active1.push_back(cumReported_Active1);  //this only applies for LR setting
      //ccumReported_Active2.push_back(cumReported_Active2);  //this only applies for LR setting
      // bbyYearReported_Passive1L.push_back(byYearReported_Passive1L); 
      // bbyYearReported_Passive2L.push_back(byYearReported_Passive2L); 
      // bbyYearReported_Passive1H.push_back(byYearReported_Passive1H); 
      // bbyYearReported_Passive2H.push_back(byYearReported_Passive2H); 
      
  
  
      // these have been stored, so reset to zero
      byYear_Incid_1L = 0;
      byYear_Incid_1H = 0;
	  
	  byYear_Diagnosis_1L = 0;
	  byYear_Diagnosis_2L = 0;
	  byYear_Diagnosis_1H = 0; 
	  byYear_Diagnosis_2H = 0;
	  
      
      byYear_Active1_fexi = 0;
      byYear_Active1_curr = 0;
      byYear_Active2_fexi = 0;
      byYear_Active2_curr = 0;
	  
	  byYear_Passive1_fexi = 0;
      byYear_Passive1_curr = 0;
      byYear_Passive2_fexi = 0;
      byYear_Passive2_curr = 0;
      
      //byYear_Passive1L_fexi = 0;
      //byYear_Passive1L_curr = 0;
      //byYear_Passive2L_fexi = 0;
      //byYear_Passive2L_curr = 0;
      //byYear_Passive1H_fexi = 0;
      //byYear_Passive1H_curr = 0;
      //byYear_Passive2H_fexi = 0;
      //byYear_Passive2H_curr = 0;
      
     // cumReported_Active1 = I1_L;
    //  cumReported_Active2 = I2_L;  //this only applies for LR setting
    // std::cout <<  I1_L << std::endl; // printing  
      cumI_1L = I1_L; //cumI_2L = 0; 
	  cumI_1H = 0; //cumI_2H = 0;
          
     // byYearReported_Passive1L = 0; 
      //  byYearReported_Passive2L = 0; 
      //  byYearReported_Passive1H = 0; 
      //  byYearReported_Passive2H = 0; 
  
      tdisc++;
      //tdisc2 = tdisc-1; //to point to period where AS is applied

      
      //   std::cout <<  r_pd1 << std::endl; // printing  
      
     // if(E_L==0 && I1_L==0 && I2_L==0 && E_H==0 && I1_H==0 && I2_H==0
     //      && EvL==0 && IvL==0 && EvH==0 && IvH==0) {break;}
    }
    
    //Determine different rates and efficiency depending on the treatment used
// For the moment, fexi is used in both AS and PD although I would maybe want to explore cases where fexi is only applied in AS,
// in such a case I will need to do some work here 'cause the boolean 'fexi' won't be enough
  
  //Determine if active screening applies or not
  //if(t<=(tdisc2+onemonth)){
  //  r_AS = rAS; // r in AS depends on omega; for the moment is the same for stage 1 and 2
  //} else {
  //  r_AS = 0;
  //}
  
  // below print to check whether annual AS correctly rate vary between years
  //if(t>(tdisc2+0.01) && t<(tdisc2+0.01001)){std::cout << t << " " << tdisc2+onemonth << " " << r_AS << std::endl;} // start new year
  //if(t>(tdisc2+0.083333333) && t<(tdisc2+0.08334)){std::cout << t << " " << tdisc2+onemonth << " " << r_AS << std::endl;} // after first month
  
  // below print to check whether annual PD rate for stg 1 correctly rate varies between years
  //if(t>(tdisc2+0.01) && t<(tdisc2+0.01001)){ // within the first month
  // std::cout << r_1L << std::endl;
  //std::cout << r_pd1 << std::endl;
  //} 
  
  // if(t>(tdisc2+0.01) && t<(tdisc2+0.01001)){std::cout << r_1L << std::endl;} // printing
  // if(t>(tdisc2+0.083333333) && t<(tdisc2+0.08334)){std::cout << r_1L << std::endl;} // printing
  
  double year_int = tdisc-tstart-1;
  double r_AS = r_ASvec[year_int];
  double r_pd1 = r_pd1vec[year_int];
  double r_pd2 = r_pd2vec[year_int];
  
  //std::cout << r_AS << std::endl;
  
  double N_L  = S_L + E_L + I1_L + I2_L;
  double N_H  = S_H + E_H + I1_H + I2_H;
  double N_aL = SaL + EaL + IaL + RaL;
  double N_aH = SaH + EaH + IaH + RaH;
  
  // Terms for biting probabilities
  double theta_vLhL = (sigma*N_L)/(sigma*(N_L+(1-xi)*N_H)+sigma_aL*N_aL);
  double theta_vLhH = (sigma*(1-xi)*N_H) / (sigma*(N_L+(1-xi)*N_H)+sigma_aL*N_aL);
  double theta_vLaL = (sigma_aL*N_aL)/(sigma*(N_L+(1-xi)*N_H)+sigma_aL*N_aL);
  double theta_vHhH = (sigma*xi*N_H)/(sigma*xi*N_H + sigma_aH*N_aH);
  double theta_vHaH = (sigma_aH*N_aH)/(sigma*xi*N_H + sigma_aH*N_aH);
  
  // Contact rates in vectors, hosts & animals
  // HUMANS
  double lambdaL  = (b*f*theta_vLhL)/N_L;
  double lambdaH1 = (b*f*theta_vLhH)/N_H; // * IvL
  double lambdaH2 = (b*f*theta_vHhH)/N_H; // * IvH
  // VECTORS
  double c1        = (f*theta_vLhL*c_h)/N_L;   // *(I1_L + I2_L)
  double c2        = (f*theta_vLhH*c_h)/N_H;   // *(I1_H + I2_H)
  double c3        = (f*theta_vLaL*c_aL)/N_aL; // *IaL
  double Lambda_vL = c1*(I1_L + I2_L)  + c2*(I1_H + I2_H) + c3*IaL;
  double c4        = (f*theta_vHhH*c_h)/N_H;   // *(I1_H + I2_H)
  double c5        = (f*theta_vHaH*c_aH)/N_aH; // *IaH
  double Lambda_vH = c4*(I1_H + I2_H)  + c5*IaH;
  // ANIMALS
  double lambda_aL = (b*f*theta_vLaL*c_aL)/N_aL;
  double lambda_aH = (b*f*theta_vHaH*c_aH)/N_aH;
  
 // initial stochastic simulation from 2012 to 2020  
  if(t<tdisc3){
   
	  if(tdisc>2016){
	  specificity = 1; //MD note: add in comment to explain
  }
	
    // Determine the total rate of each event
    //Humans Low-risk setting
    double a1 = delta*T_L;
    double a2 = mu*S_L;
    double a3 = lambdaL*IvL*S_L;
    double a4 = mu*E_L;
    double a5 = eta*E_L;
    double a6 = mu*I1_L;
    double a7 = gamma*I1_L;
    double a8_1 = r_AS*sensitivity*I1_L;
	double a8_1fp = r_AS*(1-specificity)*(S_L+E_L+T_L);  // False positivity rate 
    double a8_2 = r_pd1*I1_L;
    double a9 = mu*I2_L;
    double a10 = mu_gamma*I2_L;
    double a11_1 = r_AS*sensitivity*I2_L;
    double a11_2 = r_pd2*I2_L;
    double a12 = mu*T_L;
    double a13 = mu_t*T_L;
    
    //Humans High-risk setting
    double b1 = delta*T_H;
    double b2 = mu*S_H;
    double b3_1 = lambdaH1*IvL*S_H;
    double b3_2 = lambdaH2*IvH*S_H;
    double b4 = mu*E_H;
    double b5 = eta*E_H;
    double b6 = mu*I1_H;
    double b7 = gamma*I1_H;
    double b8 = r_pd1*I1_H;
    double b9 = mu*I2_H;
    double b10 = mu_gamma*I2_H;
    double b11 = r_pd2*I2_H;
    double b12 = mu*T_H;
    double b13 = mu_t*T_H;
    
    // Animals Low-risk setting
    double C1 = delta_a*RaL;
    double C2 = mu_aL*SaL;
    double C3 = lambda_aL*IvL*SaL;
    double C4 = mu_aL*EaL;
    double C5 = eta*EaL;
    double C6 = mu_aL*IaL;
    double C7 = gamma_aL*IaL;
    double C8 = mu_aL*RaL;
    
    // Animals High-risk setting
    double d1 = delta_a*RaH;
    double d2 = mu_aH*SaH;
    double d3 = lambda_aH*IvH*SaH;
    double d4 = mu_aH*EaH;
    double d5 = eta*EaH;
    double d6 = mu_aH*IaH;
    double d7 = gamma_aH*IaH;
    double d8 = mu_aH*RaH;
    
    // Vector Low-risk setting
    double e1 = mu_v*SvL;
    double e2 = Lambda_vL*SvL;
    double e3 = alpha*SvL;
    double e4 = mu_v*EvL;
    double e5 = nu*EvL;
    double e6 = mu_v*IvL;
    double e7 = mu_v*UvL;
    
    // Vector High-risk setting
    double f1 = mu_v*SvH;
    double f2 = Lambda_vH*SvH;
    double f3 = alpha*SvH;
    double f4 = mu_v*EvH;
    double f5 = nu*EvH;
    double f6 = mu_v*IvH;
    double f7 = mu_v*UvH;
    
    // Determine the total event rate
    double sumprobs = a1+a2+a3+a4+a5+a6+a7+a8_1+a8_1fp+a8_2+a9+a10+a11_1+a11_2+a12+a13+
      b1+b2+b3_1+b3_2+b4+b5+b6+b7+b8+b9+b10+b11+b12+b13+
      C1+C2+C3+C4+C5+C6+C7+C8+d1+d2+d3+d4+d5+d6+d7+d8+
      e1+e2+e3+e4+e5+e6+e7+f1+f2+f3+f4+f5+f6+f7;
    //std::cout << sumprobs << std::endl;  
	
    // define these constants, useful below
    double total_a = a1+a2+a3+a4+a5+a6+a7+a8_1+a8_1fp+a8_2+a9+a10+a11_1+a11_2+a12+a13;
    double total_b = b1+b2+b3_1+b3_2+b4+b5+b6+b7+b8+b9+b10+b11+b12+b13;
    double total_C = C1+C2+C3+C4+C5+C6+C7+C8;
    double total_d = d1+d2+d3+d4+d5+d6+d7+d8;
    double total_e = e1+e2+e3+e4+e5+e6+e7;
    // double total_f = f1+f2+f3+f4+f5+f6+f7;

    
    // Determine the time step at which the next event occurs + update the time
    double dt = rexp(1,sumprobs)[0];
    t += dt;
	std::cout << sumprobs << std::endl;
	std::cout << dt << std::endl;
    std::cout << t << std::endl;
    // Determine the next event to occur + update the populations
    double r = runif(1)[0];
    double p = r * sumprobs;
	std::cout << p << std::endl;
    
    double sum = a1;
    if (p <= a1) { T_L--; S_L++; sum +=a2;}
    else if (p <= a1+a2) sum += a3; // do nothing (they cancel each other)
    else if (p <= a1+a2+a3) { S_L--; E_L++; sum += a4;}
    else if (p <= a1+a2+a3+a4) { E_L--; S_L++; sum += a5;}
    else if (p <= a1+a2+a3+a4+a5) { E_L--; I1_L++; byYear_Incid_1L++; sum += a6;}
    else if (p <= a1+a2+a3+a4+a5+a6) { I1_L--; S_L++; sum += a7;}
    else if (p <= a1+a2+a3+a4+a5+a6+a7) { I1_L--; I2_L++; sum += a8_1;}
    else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1) { I1_L--; T_L++; byYear_Diagnosis_1L++; byYear_Active1_curr++; sum += a8_1fp;}
	else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1+a8_1fp) { byYear_Diagnosis_1L++; byYear_Active1_curr++; sum += a8_2;}
    else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1+a8_1fp+a8_2) { I1_L--; T_L++; byYear_Diagnosis_1L++; byYear_Passive1_curr++; sum += a9;}
    else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1+a8_1fp+a8_2+a9) { I2_L--; S_L++; sum += a10;}
    else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1+a8_1fp+a8_2+a9+a10) { I2_L--; S_L++; sum += a11_1;}
    else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1+a8_1fp+a8_2+a9+a10+a11_1) { I2_L--; T_L++; byYear_Diagnosis_2L++; byYear_Active2_curr++; sum += a11_2;}
    else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1+a8_1fp+a8_2+a9+a10+a11_1+a11_2) { I2_L--; T_L++; byYear_Diagnosis_2L++; byYear_Passive2_curr++; sum += a12;}
    else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1+a8_1fp+a8_2+a9+a10+a11_1+a11_2+a12) { T_L--; S_L++; sum += a13;}
    else if (p <= total_a) { T_L--; S_L++; sum += b1;}
    
    else if (p <= total_a+b1) { T_H--; S_H++;  sum += b2;}
    else if (p <= total_a+b1+b2) sum += b3_1;
    else if (p <= total_a+b1+b2+b3_1) { S_H--; E_H++; sum += b3_2;}
    else if (p <= total_a+b1+b2+b3_1+b3_2) { S_H--; E_H++; sum += b4;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4) { E_H--; S_H++; sum += b5;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5) { E_H--; I1_H++; byYear_Incid_1H++; sum += b6;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6) { I1_H--; S_H++; sum += b7;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7) { I1_H--; I2_H++; sum += b8;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7+b8) { I1_H--; T_H++; byYear_Diagnosis_1H++; byYear_Passive1_curr++; sum += b9;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7+b8+b9) { I2_H--; S_H++; sum += b10;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7+b8+b9+b10) { I2_H--; S_H++; sum += b11;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7+b8+b9+b10+b11) { I2_H--; T_H++; byYear_Diagnosis_2H++; byYear_Passive2_curr++; sum += b12;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7+b8+b9+b10+b11+b12) { T_H--; S_H++; sum += b13;}
    else if (p <= total_a+total_b) { T_H--; S_H++; sum += C1;}
    
    else if (p <= total_a+total_b+C1) { RaL--; SaL++; sum += C2;}
    else if (p <= total_a+total_b+C1+C2) sum += C3;
    else if (p <= total_a+total_b+C1+C2+C3) { SaL--;  EaL++; sum += C4;}
    else if (p <= total_a+total_b+C1+C2+C3+C4) { EaL--; SaL++; sum += C5;}
    else if (p <= total_a+total_b+C1+C2+C3+C4+C5) { EaL--; IaL++; sum += C6;}
    else if (p <= total_a+total_b+C1+C2+C3+C4+C5+C6) { IaL--; SaL++; sum += C7;}
    else if (p <= total_a+total_b+C1+C2+C3+C4+C5+C6+C7) { IaL--; RaL++; sum += C8;}
    else if (p <= total_a+total_b+total_C) { RaL--; SaL++; sum += d1;}
    
    else if (p <= total_a+total_b+total_C+d1) { RaH--; SaH++; sum += d2;}
    else if (p <= total_a+total_b+total_C+d1+d2) sum += d3;
    else if (p <= total_a+total_b+total_C+d1+d2+d3) { SaH--;  EaH++; sum += d4;}
    else if (p <= total_a+total_b+total_C+d1+d2+d3+d4) { EaH--; SaH++; sum += d5;}
    else if (p <= total_a+total_b+total_C+d1+d2+d3+d4+d5) { EaH--; IaH++; sum += d6;}
    else if (p <= total_a+total_b+total_C+d1+d2+d3+d4+d5+d6) { IaH--; SaH++; sum += d7;}
    else if (p <= total_a+total_b+total_C+d1+d2+d3+d4+d5+d6+d7) { IaH--; RaH++; sum += d8;}
    else if (p <= total_a+total_b+total_C+total_d) { RaH--; SaH++;  sum += e1;}
    
    else if (p <= total_a+total_b+total_C+total_d+e1) sum += e2;
    else if (p <= total_a+total_b+total_C+total_d+e1+e2) { SvL--; EvL++; sum += e3;}
    else if (p <= total_a+total_b+total_C+total_d+e1+e2+e3) { SvL--; UvL++; sum += e4;}
    else if (p <= total_a+total_b+total_C+total_d+e1+e2+e3+e4) { EvL--; SvL++; sum += e5;}
    else if (p <= total_a+total_b+total_C+total_d+e1+e2+e3+e4+e5) { EvL--; IvL++; sum += e6;}
    else if (p <= total_a+total_b+total_C+total_d+e1+e2+e3+e4+e5+e6) { IvL--; SvL++; sum += e7;}
    else if (p <= total_a+total_b+total_C+total_d+total_e) { UvL--; SvL++; sum += f1;}
    
    else if (p <= total_a+total_b+total_C+total_d+total_e+f1) sum += f2;
    else if (p <= total_a+total_b+total_C+total_d+total_e+f1+f2) { SvH--; EvH++; sum += f3;}
    else if (p <= total_a+total_b+total_C+total_d+total_e+f1+f2+f3) { SvH--; UvH++; sum += f4;}
    else if (p <= total_a+total_b+total_C+total_d+total_e+f1+f2+f3+f4) { EvH--; SvH++; sum += f5;}
    else if (p <= total_a+total_b+total_C+total_d+total_e+f1+f2+f3+f4+f5) { EvH--; IvH++; sum += f6;}
    else if (p <= total_a+total_b+total_C+total_d+total_e+f1+f2+f3+f4+f5+f6) {IvH--; SvH++;}
    else {UvH--; SvH++; }
    
    //if(t<(tdisc3-0.01) && t>(tdisc3-0.0101)){std::cout << "first section, time before tdisc3" << std::endl;} // printing
    //if(t>(tdisc3+0.01) && t<(tdisc3+0.0101)){std::cout << "first section, time after tdisc3" << std::endl;} // printing
	  
    } else{ 
  
    // Determine the total rate of each event
    //Humans Low-risk setting
    double a1 = delta*T_L;
    double a2 = mu*S_L;
    double a3 = lambdaL*IvL*S_L;
    double a4 = mu*E_L;
    double a5 = eta*E_L;
    double a6 = mu*I1_L;
    double a7 = gamma*I1_L;
	double a8_1 = (r_AS*sensitivity+r_pd1)*I1_L; // Correct diagnosis rate for stage 1
	double a8_2 = (1-psi)*(r_AS*sensitivity+r_pd1)*phi1*I1_L; //Fexi non-compliance rate - these patients remain in the I1_L class
	double a9_1 = psi*r_AS*sensitivity*phi1*I1_L; // AS using fexi to a proportion (phi1) of stage 1 cases with compliance psi
    double a9_2 = psi*r_pd1*phi1*I1_L; // PD using fexi to a proportion (phi1) of stage 1 cases, with compliance psi
    double a10_1 = r_AS*sensitivity*(1-phi1)*I1_L; // AS using the current treatment to a proportion (1-phi1) of stage 1 cases
    double a10_2 = r_pd1*(1-phi1)*I1_L; // // PD using the current treatment to a proportion (1-phi1) of stage 1 cases
    double a11 = mu*I2_L;
    double a12 = mu_gamma*I2_L;
	double a13_1 = (r_AS*sensitivity+r_pd2)*I2_L; // Diagnosis rate for stage 2
	double a13_2 = (1-psi)*(r_AS*sensitivity+r_pd2)*phi2*I2_L; //Fexi non-compliance rate - these patients remain in the I2_L class
	double a14_1 = psi*r_AS*sensitivity*phi2*I2_L; // AS using fexi to a proportion (phi2) of stage 2 cases
    double a14_2 = psi*r_pd2*phi2*I2_L; //PD using fexi to a proportion (phi2) of stage 2 cases, with efficacy e2 (and compliance omega2, hidden in r_pd2Fexi)
    double a15_1 = r_AS*sensitivity*(1-phi2)*I2_L; // AS using the current treatment to a proportion (1-phi2) of stage 2 cases
    double a15_2 = r_pd2*(1-phi2)*I2_L; // PD using the current treatment to a proportion (1-phi2) of stage 2 cases
    double a16 = mu*T_L;
    double a17 = mu_t*T_L;
    
    
    //Humans High-risk setting
    double b1 = delta*T_H;
    double b2 = mu*S_H;
    double b3_1 = lambdaH1*IvL*S_H;
    double b3_2 = lambdaH2*IvH*S_H;
    double b4 = mu*E_H;
    double b5 = eta*E_H;
    double b6 = mu*I1_H;
    double b7 = gamma*I1_H;
	double b8_1 = r_pd1*I1_H;
	double b8_2 = (1-psi)*r_pd1*phi1*I1_H;
    double b9 = psi*r_pd1*phi1*I1_H; // PD using fexi to a proportion (phi1) of stage 1 cases with compliance psi
    double b10 = r_pd1*(1-phi1)*I1_H; // PD using the current treatment to a proportion (1-phi2) of stage 2 cases with compliance psi
    double b11 =  mu*I2_H;
    double b12 = mu_gamma*I2_H;
	double b13_1 = r_pd2*I2_H;
	double b13_2 = (1-psi)*r_pd2*phi2*I2_H;
    double b14 = psi*r_pd2*phi2*I2_H; //PD using fexi to a proportion (phi2) of stage 2 cases, with efficacy e2 (and compliance omega2, hidden in r_pd2Fexi)
    double b15 = r_pd2*(1-phi2)*I2_H; // PD using the current treatment to a proportion (1-phi2) of stage 2 cases    
    double b16 = mu*T_H;
    double b17 = mu_t*T_H;

    // Animals Low-risk setting
    double C1 = delta_a*RaL;
    double C2 = mu_aL*SaL;
    double C3 = lambda_aL*IvL*SaL;
    double C4 = mu_aL*EaL;
    double C5 = eta*EaL;
    double C6 = mu_aL*IaL;
    double C7 = gamma_aL*IaL;
    double C8 = mu_aL*RaL;

    // Animals High-risk setting
    double d1 = delta_a*RaH;
    double d2 = mu_aH*SaH;
    double d3 = lambda_aH*IvH*SaH;
    double d4 = mu_aH*EaH;
    double d5 = eta*EaH;
    double d6 = mu_aH*IaH;
    double d7 = gamma_aH*IaH;
    double d8 = mu_aH*RaH;

    // Vector Low-risk setting
    double e1 = mu_v*SvL;
    double e2 = Lambda_vL*SvL;
    double e3 = alpha*SvL;
    double e4 = mu_v*EvL;
    double e5 = nu*EvL;
    double e6 = mu_v*IvL;
    double e7 = mu_v*UvL;

    // Vector High-risk setting
    double f1 = mu_v*SvH;
    double f2 = Lambda_vH*SvH;
    double f3 = alpha*SvH;
    double f4 = mu_v*EvH;
    double f5 = nu*EvH;
    double f6 = mu_v*IvH;
    double f7 = mu_v*UvH;
    
    // define these constants, useful below
    double total_a = a1+a2+a3+a4+a5+a6+a7+a8_2+a9_1+a9_2+a10_1+a10_2+a11+a12+a13_2+a14_1+a14_2+a15_1+a15_2+a16+a17; //a8_1 and a13_1 not included as they are the diagnosis rates
    double total_b = b1+b2+b3_1+b3_2+b4+b5+b6+b7+b8_2+b9+b10+b11+b12+b13_2+b14+b15+b16+b17; ////b8_1 and b13_1 not included as they are the diagnosis rates
    double total_C = C1+C2+C3+C4+C5+C6+C7+C8;
    double total_d = d1+d2+d3+d4+d5+d6+d7+d8;
    double total_e = e1+e2+e3+e4+e5+e6+e7;
    // double total_f = f1+f2+f3+f4+f5+f6+f7;

    // Determine the total event rate
    double sumprobs = total_a + total_b + total_C + total_d + total_e + 
      f1+f2+f3+f4+f5+f6+f7;
    
    // Determine the time step at which the next event occurs + update the time
    double dt = rexp(1,sumprobs)[0];
    t += dt;

    // Determine the next event to occur + update the populations
    double r = runif(1)[0];
    double p = r * sumprobs;
    
    double sum = a1;
    if (p <= a1) { T_L--; S_L++; sum +=a2;}
    else if (p <= a1+a2) sum += a3; // do nothing (they cancel each other)
    else if (p <= a1+a2+a3) { S_L--; E_L++; sum += a4;}
    else if (p <= a1+a2+a3+a4) { E_L--; S_L++; sum += a5;}
    else if (p <= a1+a2+a3+a4+a5) { E_L--; I1_L++; byYear_Incid_1L++; cumI_1L++; sum += a6;}
    else if (p <= a1+a2+a3+a4+a5+a6) { I1_L--; S_L++; sum += a7;}
    else if (p <= a1+a2+a3+a4+a5+a6+a7) { I1_L--; I2_L++; sum += a8_1;}
	else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1) { byYear_Diagnosis_1L++; //Diagnosed, now need to classify by treatment type
		if (p <= a1+a2+a3+a4+a5+a6+a7+a9_1) { I1_L--; T_L++; byYear_Active1_fexi++;}
		else if (p <= a1+a2+a3+a4+a5+a6+a7+a9_1+a9_2) { I1_L--; T_L++; byYear_Passive1_fexi++;}
		else if (p <= a1+a2+a3+a4+a5+a6+a7+a9_1+a9_2+a10_1)	{ I1_L--; T_L++; byYear_Active1_curr++;}
		else if (p <= a1+a2+a3+a4+a5+a6+a7+a9_1+a9_2+a10_1+a10_2) { I1_L--; T_L++; byYear_Passive1_curr++;}
		sum += a8_2 + a9_1 + a9_2 + a10_1 + a10_2;} //else do nothing, diagnosed but not treated
    else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1+a11) { I2_L--; S_L++; sum += a12;}
    else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1+a11+a12) { I2_L--; S_L++; sum += a13_1;}
	else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1+a11+a12+a13_1) { byYear_Diagnosis_2L++; 
		if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1+a11+a12+a14_1) { I2_L--; T_L++; byYear_Active2_fexi++;}
		else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1+a11+a12+a14_1+a14_2) { I2_L--; T_L++; byYear_Passive2_fexi++;}
		else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1+a11+a12+a14_1+a14_2+a15_1) { I2_L--; T_L++; byYear_Active2_curr++;} 
		else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1+a11+a12+a14_1+a14_2+a15_1+a15_2) { I2_L--; T_L++; byYear_Passive2_curr++;}
		sum += a13_2 + a14_1 + a14_2 + a15_1 + a15_2;}
    else if (p <= a1+a2+a3+a4+a5+a6+a7+a8_1+a11+a12+a13_1+a16) { T_L--; S_L++; sum += a17;}
    else if (p <= total_a) { T_L--; S_L++; sum += b1;}
    
    else if (p <= total_a+b1) { T_H--; S_H++;  sum += b2;}
    else if (p <= total_a+b1+b2) sum += b3_1;
    else if (p <= total_a+b1+b2+b3_1) { S_H--; E_H++; sum += b3_2;}
    else if (p <= total_a+b1+b2+b3_1+b3_2) { S_H--; E_H++; sum += b4;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4) { E_H--; S_H++; sum += b5;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5) { E_H--; I1_H++; byYear_Incid_1H++; cumI_1H++; sum += b6;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6) { I1_H--; S_H++; sum += b7;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7) { I1_H--; I2_H++; sum += b8_1;}
	else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7+b8_1) { byYear_Diagnosis_1H++; 
		if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7+b9) { I1_H--; T_H++; byYear_Passive1_fexi++;}
		else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7+b9+b10) { I1_H--; T_H++; byYear_Passive1_curr++;}
		sum += b8_2 + b9 + b10;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7+b8_1+b11) { I2_H--; S_H++; sum += b12;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7+b8_1+b11+b12) { I2_H--; S_H++; sum += b13_1;}
	else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7+b8_1+b11+b12+b13_1) { byYear_Diagnosis_2H++; 
		if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7+b8_1+b11+b12+b14) { I2_H--; T_H++; byYear_Passive2_fexi++;}
		else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7+b8_1+b11+b12+b14+b15) { I2_H--; T_H++; byYear_Passive2_curr++;}
		sum += b13_2 + b14 + b15;}
    else if (p <= total_a+b1+b2+b3_1+b3_2+b4+b5+b6+b7+b8_1+b11+b12+b13_1+b16) { T_H--; S_H++; sum += b17;}
    else if (p <= total_a+total_b) { T_H--; S_H++; sum += C1;}

    else if (p <= total_a+total_b+C1) { RaL--; SaL++; sum += C2;}
    else if (p <= total_a+total_b+C1+C2) sum += C3;
    else if (p <= total_a+total_b+C1+C2+C3) { SaL--;  EaL++; sum += C4;}
    else if (p <= total_a+total_b+C1+C2+C3+C4) { EaL--; SaL++; sum += C5;}
    else if (p <= total_a+total_b+C1+C2+C3+C4+C5) { EaL--; IaL++; sum += C6;}
    else if (p <= total_a+total_b+C1+C2+C3+C4+C5+C6) { IaL--; SaL++; sum += C7;}
    else if (p <= total_a+total_b+C1+C2+C3+C4+C5+C6+C7) { IaL--; RaL++; sum += C8;}
    else if (p <= total_a+total_b+total_C) { RaL--; SaL++; sum += d1;}

    else if (p <= total_a+total_b+total_C+d1) { RaH--; SaH++; sum += d2;}
    else if (p <= total_a+total_b+total_C+d1+d2) sum += d3;
    else if (p <= total_a+total_b+total_C+d1+d2+d3) { SaH--;  EaH++; sum += d4;}
    else if (p <= total_a+total_b+total_C+d1+d2+d3+d4) { EaH--; SaH++; sum += d5;}
    else if (p <= total_a+total_b+total_C+d1+d2+d3+d4+d5) { EaH--; IaH++; sum += d6;}
    else if (p <= total_a+total_b+total_C+d1+d2+d3+d4+d5+d6) { IaH--; SaH++; sum += d7;}
    else if (p <= total_a+total_b+total_C+d1+d2+d3+d4+d5+d6+d7) { IaH--; RaH++; sum += d8;}
    else if (p <= total_a+total_b+total_C+total_d) { RaH--; SaH++;  sum += e1;}

    else if (p <= total_a+total_b+total_C+total_d+e1) sum += e2;
    else if (p <= total_a+total_b+total_C+total_d+e1+e2) { SvL--; EvL++; sum += e3;}
    else if (p <= total_a+total_b+total_C+total_d+e1+e2+e3) { SvL--; UvL++; sum += e4;}
    else if (p <= total_a+total_b+total_C+total_d+e1+e2+e3+e4) { EvL--; SvL++; sum += e5;}
    else if (p <= total_a+total_b+total_C+total_d+e1+e2+e3+e4+e5) { EvL--; IvL++; sum += e6;}
    else if (p <= total_a+total_b+total_C+total_d+e1+e2+e3+e4+e5+e6) { IvL--; SvL++; sum += e7;}
    else if (p <= total_a+total_b+total_C+total_d+total_e) { UvL--; SvL++; sum += f1;}

    else if (p <= total_a+total_b+total_C+total_d+total_e+f1) sum += f2;
    else if (p <= total_a+total_b+total_C+total_d+total_e+f1+f2) { SvH--; EvH++; sum += f3;}
    else if (p <= total_a+total_b+total_C+total_d+total_e+f1+f2+f3) { SvH--; UvH++; sum += f4;}
    else if (p <= total_a+total_b+total_C+total_d+total_e+f1+f2+f3+f4) { EvH--; SvH++; sum += f5;}
    else if (p <= total_a+total_b+total_C+total_d+total_e+f1+f2+f3+f4+f5) { EvH--; IvH++; sum += f6;}
    else if (p <= total_a+total_b+total_C+total_d+total_e+f1+f2+f3+f4+f5+f6) {IvH--; SvH++;}
    else {UvH--; SvH++; }

    //  if (t<=(tdisc+tole)) {cumI_1L+=I1_L;}
    // std::cout << cumI_1L << std::endl; // printing
    
    //cumI_2L+=I2_L; cumI_1H+=I1_H; cumI_2H+=I2_H; //these are annual cumulative values
          
          // aqui tratando de hacer cumulative para el primer mes, asi se cuanto hay reported via Active screening          
          // if(t<=tdisc2+onemonth) { 
          //       cumReported_Active1+ =I1_L; cumReported_Active2+=I2_L; 
          // }  
          //if(t<(tdisc3-0.01) && t>(tdisc3-0.0101)){std::cout << "second section, time before tdisc3" << std::endl;} // printing
          //if(t>(tdisc3+0.01) && t<(tdisc3+0.0101)){std::cout << "second section, time after tdisc3" << std::endl;} // printing

  }        
  } while (t<=(maxyears+2000+tole));// && (I1_L>0 || I1_H>0));
  // 
  // std::cout << t << std::endl; // printing 
      
      
  //create a list of (state) vectors and plug into a data frame structure
  int ncols =15; 
  
  Rcpp::List long_list(ncols);
  
  long_list[0] = tt;
  /* //humans
  long_list[1] = SS_L;
  long_list[2] = EE_L;
  long_list[3] = II1_L;
  long_list[4] = II2_L;
  long_list[5] = TT_L;
  
  long_list[6] = SS_H;
  long_list[7] = EE_H;
  long_list[8] = II1_H;
  long_list[9] = II2_H;
  long_list[10] = TT_H;
  
  //vectors
  long_list[11] = SSvL;
  long_list[12] = EEvL;
  long_list[13] = IIvL;
  long_list[14] = UUvL;
  
  long_list[15] = SSvH;
  long_list[16] = EEvH;
  long_list[17] = IIvH;
  long_list[18] = UUvH;
  
  //animals
  long_list[19] = SSaL;
  long_list[20] = EEaL;
  long_list[21] = IIaL;
  long_list[22] = RRaL;
  
  long_list[23] = SSaH;
  long_list[24] = EEaH;
  long_list[25] = IIaH;
  long_list[26] = RRaH; */
  
  long_list[1] = bbyYear_Incid_1L;
  long_list[2] = bbyYear_Incid_1H;
  long_list[3] = bbyYear_Diagnosis_1L;
  long_list[4] = bbyYear_Diagnosis_2L;
  long_list[5] = bbyYear_Diagnosis_1H;
  long_list[6] = bbyYear_Diagnosis_2H;
  long_list[7] = bbyYear_Active1_fexi;
  long_list[8] = bbyYear_Active1_curr;
  long_list[9] = bbyYear_Active2_fexi;
  long_list[10] = bbyYear_Active2_curr;
  long_list[11] = bbyYear_Passive1_fexi;
  long_list[12] = bbyYear_Passive1_curr;
  long_list[13] = bbyYear_Passive2_fexi;
  long_list[14] = bbyYear_Passive2_curr;  
  
  //long_list[11] = bbyYear_Passive1L_fexi;
  //long_list[12] = bbyYear_Passive1L_curr;
  //long_list[13] = bbyYear_Passive2L_fexi;
  //long_list[14] = bbyYear_Passive2L_curr;
  //long_list[15] = bbyYear_Passive1H_fexi;
  //long_list[16] = bbyYear_Passive1H_curr;
  //long_list[17] = bbyYear_Passive2H_fexi;
  //long_list[18] = bbyYear_Passive2H_curr;

   // std::vector<char> namevec(tt.size());
  // char colNames = params["colNames"]; 
  //std::string fname = Rcpp::as<std::string>(colNames); 
  // Add colnames
  //long_list.attr("names") = Rcpp::as<std::string>(colNames);// colNames;
  
  // Coerce list to data.frame
  long_list.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, tt.size());
  long_list.attr("class") = "data.frame";
  return long_list;
};


