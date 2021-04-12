#define TMB_LIB_INIT R_init_mpmm
#include <TMB.hpp>
#include "TMB/sub/mp.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(model_name);
  if (model_name == "mp") {
    return mp(this);
  } else {
    error ("Unknown model_name");
  }
  return 0;
}
