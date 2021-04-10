#include <cmath>
#ifndef mp_hpp
#define mp_hpp 1

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace density;

// modified from glmmTMB
// need to link to these via glmmTMB pkg but can't currently because they do not export C++
template <class Type>
struct per_term_info {
  // Input from R
  int blockCode;     // Code that defines structure
  int blockSize;     // Size of one block
  int blockReps;     // Repeat block number of times
  int blockNumTheta; // Parameter count per block
  matrix<Type> dist;
  // Report output
  matrix<Type> corr;
  vector<Type> sd;
};
// need to link to these via glmmTMB pkg
template <class Type>
struct terms_t : vector<per_term_info<Type> > {
  terms_t(SEXP x){
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP y = VECTOR_ELT(x, i);    // y = x[[i]]
      int blockCode = (int) REAL(getListElement(y, "blockCode", &isNumericScalar))[0];
      int blockSize = (int) REAL(getListElement(y, "blockSize", &isNumericScalar))[0];
      int blockReps = (int) REAL(getListElement(y, "blockReps", &isNumericScalar))[0];
      int blockNumTheta = (int) REAL(getListElement(y, "blockNumTheta", &isNumericScalar))[0];
      (*this)(i).blockCode = blockCode;
      (*this)(i).blockSize = blockSize;
      (*this)(i).blockReps = blockReps;
      (*this)(i).blockNumTheta = blockNumTheta;
    }
  }
};

enum valid_covStruct {
  diag_covstruct = 0,
  us_covstruct   = 1
};
// would be better to have AR(1) corrStruct here to deal with serial dependence in animal tracks
// need to link to these via glmmTMB pkg but can't currently because they do not export C++
template <class Type>
Type termwise_nll(array<Type> &U, vector<Type> theta, per_term_info<Type>& term) {
  Type ans = 0;
  if (term.blockCode == diag_covstruct){
    // case: diag_covstruct
    vector<Type> sd = exp(theta);
    for(int i = 0; i < term.blockReps; i++){
      ans -= dnorm(vector<Type>(U.col(i)), Type(0), sd, true).sum();
    }
    term.sd = sd;
  }
  else if (term.blockCode == us_covstruct){
    // case: us_covstruct
    int n = term.blockSize;
    vector<Type> logsd = theta.head(n);
    vector<Type> corr_transf = theta.tail(theta.size() - n);
    vector<Type> sd = exp(logsd);
    density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
    density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
    }
    term.corr = nldens.cov();
    term.sd = sd;
  }
  else error("covStruct not implemented!");
  return ans;
}
// need to link to these via glmmTMB pkg but can't currently because they do not export C++
template <class Type>
Type allterms_nll(vector<Type> &u, vector<Type> theta,
                  vector<per_term_info<Type> >& terms) {
  Type ans = 0;
  int upointer = 0;
  int tpointer = 0;
  int nr, np = 0, offset;
  for(int i=0; i < terms.size(); i++){
    nr = terms(i).blockSize * terms(i).blockReps;
    // Note: 'blockNumTheta=0' ==> Same parameters as previous term.
    bool emptyTheta = ( terms(i).blockNumTheta == 0 );
    offset = ( emptyTheta ? -np : 0 );
    np     = ( emptyTheta ?  np : terms(i).blockNumTheta );
    vector<int> dim(2);
    dim << terms(i).blockSize, terms(i).blockReps;
    array<Type> useg( &u(upointer), dim);
    vector<Type> tseg = theta.segment(tpointer + offset, np);
    ans += termwise_nll(useg, tseg, terms(i));
    upointer += nr;
    tpointer += terms(i).blockNumTheta;
  }
  return ans;
}

// mpmm model code
template <class Type>
Type mp(objective_function<Type>* obj) {

  DATA_MATRIX(X);                   // fixed effects design matrix
  DATA_SPARSE_MATRIX(Z);            // random effects design sparse matrix
  DATA_ARRAY(ll);                  // locations
  DATA_FACTOR(idx);                 // cumsum of number of locations for each animal
  DATA_VECTOR(di);                  // di is time interval between ll_i and ll_{i-1}
  DATA_INTEGER(A);                // number of animals
  // Define covariance structure for the conditional model
  DATA_STRUCT(terms, terms_t);

  // For one-step-ahead resisuals
  DATA_ARRAY_INDICATOR(keep, ll);

  PARAMETER_VECTOR(lg);		          // move autocorrelation parameter (gamma on link scale)
  PARAMETER_VECTOR(beta);           // fixed regression coefficients (link scale)
  PARAMETER_VECTOR(b);              // random intercept & slope terms
  PARAMETER_VECTOR(l_sigma);	    // innovation variance (log scale)
  PARAMETER(l_rho);                  // innovation correlation (logit scale)
  PARAMETER(l_sigma_g);           // logistic scale parameter of rw on lg (log scale)
  PARAMETER_VECTOR(theta);          // covariance parameters

// Backtransform parameters from link scale
vector<Type> gamma = Type(1.0) / (Type(1.0) + exp(-lg));
vector<Type> sigma = exp(l_sigma);
Type rho = Type(2.0) / (Type(1.0) + exp(-l_rho)) - Type(1.0);
Type sigma_g = exp(l_sigma_g);

// 2x2 covariance matrix for innovations - diag in loop
matrix<Type> cov(2, 2);

Type jnll = 0.0;
vector<Type> mu(2);

int i,j;

// Random intercept & slope(s)
jnll += allterms_nll(b, theta, terms);

// linear predictor
vector<Type> eta = X * beta + Z * b;

  for(i = 0; i < A; ++i) {
    jnll -= dnorm(lg(idx(i)), eta(idx(i)), sigma_g, true);
    for(j = (idx(i)+1); j < idx(i+1); ++j) {
      jnll -= dnorm(lg(j), eta(j), di(j) * sigma_g, true);
    }

    for(j = (idx(i)+2); j < idx(i+1); ++j){
      mu = ll.matrix().row(j) - ll.matrix().row(j-1) - gamma(j) * (di(j)/di(j-1)) * (ll.matrix().row(j-1) - ll.matrix().row(j-2));  // first diff RW on locations
      // var-cov depends on time interval
      cov.setZero();
      cov(0,0) = sigma(0) * sigma(0) * di(j) * di(j);
      cov(1,1) = sigma(1) * sigma(1) * di(j) * di(j);
      cov(1,0) = sigma(0) * sigma(1) * di(j) * rho;
      cov(0,1) = cov(1,0);
      MVNORM_t<Type> nll_dens(cov);
      jnll += nll_dens(mu, keep.matrix().row(j));
    }
  }

  // Report variance components
  vector<matrix<Type> > corr(terms.size());
  vector<vector<Type> > sd(terms.size());
  for(int i=0; i<terms.size(); i++){
    if(terms(i).blockNumTheta > 0){
      corr(i) = terms(i).corr;
      sd(i) = terms(i).sd;
    }
  }

  REPORT(corr);
  REPORT(sd);
  ADREPORT(sigma_g);
  ADREPORT(sigma);
  ADREPORT(rho);
  ADREPORT(beta);

  return jnll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
