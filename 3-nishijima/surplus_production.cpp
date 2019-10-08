// glmm with Poisson distribution 

#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA //
  DATA_VECTOR(Catch);
  DATA_VECTOR(cpue);
  DATA_INTEGER(SP_type); //0:Scheffer, 1:Fox
  
  // PARAMETER //
  PARAMETER(log_r); // fixed effect
  PARAMETER(log_K);
  PARAMETER(log_sigma_pro);
  PARAMETER(log_q); //catchability
  PARAMETER(log_sigma_obs);
  PARAMETER_VECTOR(log_B);
  
  // Parameter transformation
  Type r = exp(log_r);
  Type K = exp(log_K);
  Type sigma_pro = exp(log_sigma_pro);
  // Type q = exp(log_q);
  Type sigma_obs = exp(log_sigma_obs);
  vector<Type> B = exp(log_B);
  
  Type nll=0;
  
  for(int i=1;i<Catch.size();i++){ //Process likelihood
    Type pred_B = B(i-1);
    if(SP_type==0){ //Scheffer
      pred_B += r*B(i-1)*(1-B(i-1)/K)-Catch(i-1);
    }else{ // Fox
      pred_B += r*B(i-1)*(1-log_B(i-1)/log_K)-Catch(i-1);
    }
    nll -= dnorm(log_B(i),log(pred_B),sigma_pro,true);
  }
  
  // observation likelihood
  nll -= sum(dnorm(log(cpue),log_q+log_B,sigma_obs,true));
  
  ADREPORT(B);
  return nll;
}