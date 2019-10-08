// surplus production model

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
    Type pred_SP = 0;
    if(SP_type==0){ //Scheffer
      pred_SP += r*B(i-1)*(1-B(i-1)/K);
    }else{ // Fox
      pred_SP += r*B(i-1)*(log_K-log_B(i-1));
    }
    nll -= dnorm(log(B(i)-B(i-1)+Catch(i-1)),log(pred_SP),sigma_pro,true);
  }
  
  // observation likelihood
  nll -= sum(dnorm(log(cpue),log_q+log_B,sigma_obs,true));
  
  ADREPORT(B);
  return nll;
}