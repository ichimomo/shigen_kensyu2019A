// glmm with Poisson distribution 

#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA //
  DATA_VECTOR(Y);
  DATA_VECTOR(X);
  
  // PARAMETER //
  PARAMETER_VECTOR(beta); // fixed effect
  PARAMETER_VECTOR(epsilon); // random effect
  
  vector<Type> r(X.size());

  Type nll=0;
  for(int i=0;i<X.size();i++){
    r(i) = beta(CppAD::Integer(Type(0)))+beta(CppAD::Integer(Type(1)))*X(i)+epsilon(i);
    nll -= dpois(Y(i),exp(r(i)),true); // poisson distribution
    nll -= dnorm(epsilon(i),Type(0),beta(CppAD::Integer(Type(2))),true); // random effect
  }
  return nll;
}