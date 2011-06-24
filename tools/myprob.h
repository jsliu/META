
#include "./newmat/newmatio.h"
#include "./newmat/newmatap.h"
#include "./newran/newran.h"
#include "./cprob/cprob.h"


void rdirichlet(double *ans, double *probs, int *l);
double ddirichlet(double *x, double *pars, int *l);
double dnorm(double x, double *mu, double *sd);
double genchi(int df);
int rmultinom_unif(int n);
double genbet(float aa, float bb);
double rgamma(double shape, double rate);
double dgamma(double x, double a, double b, int flag);
void matrix_exp(vector<vector<double> > *mat, double m);
void samp(int N, int M, int p, vector<int> *vec);
void samp1(int N, int M, int p, vector<int> *vec);
int countLines (istream& in);
double log_inverse_Wishart_density(SymmetricMatrix *W, double v, SymmetricMatrix *S);
double log_inverse_Wishart_density_2D(vector<double> *W, double v, vector<double> *S);
double log_mvt_norm_density(ColumnVector *x, ColumnVector *mu, SymmetricMatrix *V);
double covvec(vector<double> *vv1, vector<double> *vv2, double mu1, double mu2, int n);
double covvec1(vector<double> *vv1, vector<double> *vv2, int n);
double varvec1(vector<double> *vv, int n);
double varvec(vector<double> *vv, double mu, int n);
double meanvec(vector<double> *vv, int n);
void rwish(SymmetricMatrix *Ans, int v, SymmetricMatrix *S, bool inv);
double log_dbeta(double x, double a, double b);
double log_dgamma(double x, double a, double b);
double log_dnorm(double x, double mu, double var);
double log_dscaled_Inv_Chi_Squared(double x, double v, double s2);
int rmultinom(vector<double> *probs);
map<double, vector<int>, greater<double> >::iterator rmultinom(map<double, vector<int>, greater<double> > *m1, double total);
void sequential_mean(double*, double*, double);
void sequential_variance(double*, double*, double*, double);
void sequential_sd(double*, double*, double*, double);
void sequential_mean_vec(double*, double*, double, int);
void sequential_variance_vec(double*, double*, double*, double, int);
void sequential_sd_vec(double*, double*, double*, double, int);



void rdirichlet(double *ans, double *probs, int *l){
	
	/* generate a single draw from a Dirichlet(probs) distribution */
	int i;
	double t = 0.0;
	
	for(i = 0; i < *l; i++) {
		*(ans + i) = rgamma(*(probs + i), 1.0);
		t += *(ans + i);
	}
	for(i = 0; i < *l; i++) *(ans + i) = (*(ans + i) + (0.01 * t)) / (t * (1 + (*l * 0.01)));
	
}

double ddirichlet(double *x, double *pars, int *l) {

	/* returns the log-density of x from a Dirichlet distribution with parameter vector pars */
	int i;
	double t1 = 0.0, log_dens = 0.0;
		
	for(i = 0; i < *l; i++) {
		t1 += *(pars + i);
		log_dens += (*(pars + i) - 1) * log(*(x + i)) - lgamma(*(pars + i));
	}
	
	log_dens += lgamma(t1);
			
	return(log_dens);
}


double dnorm(double x, double *mu, double *sd) {

	double a;

	a = exp(-0.5 * (x - *mu) * (x - *mu) / (*sd * *sd)) / (sqrt(2 * PI) * *sd);

	return a;
}



double log_dscaled_Inv_Chi_Squared(double x, double v, double s2) { // Gelman, Carlin, Rubin and Stone definition
  
  
  double l = 0.0;
  
  l = 0.5 * v * log(s2 * v / 2) - lgamma(v / 2) - (v / 2 + 1) * log(x) - (v * s2 / (2 * x));
  
  return(l);
}

double log_dbeta(double x, double a, double b) {
  
  
  double l = 0.0;
  
  l = (a - 1) * log(x) + (b - 1) * log(1 - x) + lgamma(a + b) - lgamma(a) - lgamma(b);
  
  return(l);
}


double log_dgamma(double x, double a, double b) {
	  
  double l = a * log(b) - lgamma(a) + (a - 1) * log(x) - b * x;

  return l;

}


double log_dnorm(double x, double mu, double var) {

  double l = -0.5 * ((x - mu) * (x - mu) / var + log(2 * PI * var));
	  
  return l;

}

void rwish(SymmetricMatrix *Ans, int v, SymmetricMatrix *S, bool inv) {


  int k, l, p = (*S).Nrows();
  UpperTriangularMatrix L, Z;
  SymmetricMatrix S1;
  Matrix X;
  Normal snorm;
 
  /*   for(k = 1; k <= p; k++) { */
  /*        for(l = 1; l <= p; l++) */
  /* 	 cout << (*S)(k,l) << " "; */
  /*        cout << endl; */
  /*      } */
  /*      cout << endl; */

  if(inv) {
    S1 = (*S).i();
    L = (Cholesky(S1)).t();
  } else {
    L = (Cholesky(*S)).t();
  }

  /*  for(k = 1; k <= p; k++) { */
  /*        for(l = k; l <= p; l++) */
  /* 	 cout << L(k,l) << " "; */
  /*        cout << endl; */
  /*      } */
  /*      cout << endl; */

  Z = L;
  for(k =1; k <= p; k++) {
    Z(k,k) = sqrt(genchi(v - k + 1));
    for(l =k+1; l <= p; l++) Z(k,l) = snorm.Next();
  }
  /*  for(k = 1; k <= p; k++) { */
  /*        for(l = k; l <= p; l++) */
  /* 	 cout << Z(k,l) << " "; */
  /*        cout << endl; */
  /*      } */
  /*      cout << endl; */

  L = Z * L;
  /*  for(k = 1; k <= p; k++) { */
  /*        for(l = k; l <= p; l++) */
  /* 	 cout << L(k,l) << " "; */
  /*        cout << endl; */
  /*      } */
  /*      cout << endl; */

  X = L.t() * L;
  if(inv) X= X.i();

  *Ans << X;
     
  /*  for(k = 1; k <= p; k++) { */
  /*        for(l = 1; l <= p; l++)  */
  /* 	 cout << (*Ans)(k,l) << " "; */
  /*        cout << endl; */
  /*      } */
  /*      cout << endl; */
}




double meanvec(vector<double> *vv, int n) {
	  
  double a = 0.0;
             
  for(int i = 0; i < n; i++) a += (*vv)[i];
  return a / n;
	  
}
	
double varvec(vector<double> *vv, double mu, int n) {
	  
  double a = 0.0;
           
  for(int i = 0; i < n; i++) a += ((*vv)[i] - mu) * ((*vv)[i] - mu);
  return a / (n - 1);
	  
}

double varvec1(vector<double> *vv, int n) {
	  
  double a = 0.0, mu = meanvec(vv, n);
               
  for(int i = 0; i < n; i++) a += ((*vv)[i] - mu) * ((*vv)[i] - mu);
  return a / (n - 1);
	  
}
	
double covvec1(vector<double> *vv1, vector<double> *vv2, int n) {
	  
  double a = 0.0;
  double mu1 = meanvec(vv1, n);
  double mu2 = meanvec(vv2, n);

  for(int i = 0; i < n; i++) a += ((*vv1)[i] - mu1) * ((*vv2)[i] - mu2);
  return a / (n - 1);
	  
}
	
double covvec(vector<double> *vv1, vector<double> *vv2, double mu1, double mu2, int n) {
	  
  double a = 0.0;
  for(int i = 0; i < n; i++) a += ((*vv1)[i] - mu1) * ((*vv2)[i] - mu2);
  return a / (n - 1);
	  
}
	






	
double log_mvt_norm_density(ColumnVector *x, ColumnVector *mu, SymmetricMatrix *V) {
	  
  double a;
  a = -0.5 * (log( (*V).Determinant() ) - (*mu).ncols() * log(2.0 * PI) - ((((*x - *mu).t()) * (*V).i()) * (*x - *mu)).AsScalar());
  return a;
	  
}
double log_inverse_Wishart_density_2D(vector<double> *W, double v, vector<double> *S) {
	  
  int i, j, k;
  double a, s, w;
  a = 0.0;

  /*  if(print) { */
  /* 	    cout << "log_inverse_Wishart_density_2D" << endl; */
  /* 	    cout << v << endl; */
  /* 	    for(j = 0; j < 3; j++) cout << (*W)[j] << " "; */
  /* 	    cout << endl; */
  /* 	    for(j = 0; j < 3; j++) cout << (*S)[j] << " "; */
  /* 	    cout << endl; */
  /* 	  } */
  for(j = 0; j < 2; j++) a -= lgamma(0.5 * (v -j));
  a -= ((0.5 * v * 2) * log(2.0) + (0.25 * 2 * (2 - 1)) * log(PI));
  s = (*S)[0] * (*S)[2] - (*S)[1] * (*S)[1]; 
  w = (*W)[0] * (*W)[2] - (*W)[1] * (*W)[1]; 
  a += (0.5 * v) * log(s);
  a -= (0.5 * (v + 2 + 1)) * log(w);
  a -= 0.5 * (((*S)[0] * (*W)[2] + (*S)[2] * (*W)[0] - 2.0 *(*S)[1] * (*W)[1]) / w);
  return a;
	  
}
	
double log_inverse_Wishart_density(SymmetricMatrix *W, double v, SymmetricMatrix *S) {
	  
  int i, j, k;
  double a, s, w;
  SymmetricMatrix s1;
  a = 0.0;
  k = (*W).nrows();
  for(j = 0; j < k; j++) a -= lgamma(0.5 * (v -j));
  a -= ((0.5 * v * k) * log(2.0) + (0.25 * k * (k - 1)) * log(PI));
  s = ((*S).log_determinant()).log_value();
  w = ((*W).log_determinant()).log_value();
  a += (0.5 * v) * s;
  a -= (0.5 * (v + k + 1)) * w;
  s1 = (*S) * ((*W).i());
  a -= 0.5 * s1.trace();
  return a;
	  
}
	
void samp(int N, int M, int p, vector<int> *vec) { 

  // random sampling of p integers without replacement from 0:(N-1) with first M numbers set to 0:(M-1)
  int i, j, tmp;
  
  if(p < M || p > N) cout << "need p >= M and p <= N in samp function" << endl;
  
  set<int> pool;
  set<int>::iterator pos;
  
  for(i = 0; i < M; i++) {
    (*vec)[i] = i;
  }
  
  for(i = M; i < N; i++) 
    pool.insert(i);

  for(i = 0; i < (p - M); i++) {
    
    tmp = rmultinom_unif(N - M - i) - 1;
    pos = pool.begin();
    advance(pos, tmp);
    (*vec)[M+i] = *pos;
    pool.erase(pos);
  }

}

void samp1(int N, int M, int p, vector<int> *vec) { 
  
  // random sampling of p integers without replacement from 0:(N-1) with first M numbers set to 0:(M-1) and numbers occuring in pairs
  int i, j, tmp;
  
  if(p < M || p > N) cout << "need p >= M and p <= N in samp function" << endl;
  
  set<int> pool;
  set<int>::iterator pos;
  
  for(i = 0; i < M; i++) 
    (*vec)[i] = i;
  
  for(i = M; i < N; i+=2) 
    pool.insert(i);

  for(i = 0; i < (N - M) / 2; i++) {
    
    tmp = rmultinom_unif(pool.size()) - 1;
    pos = pool.begin();
    advance(pos, tmp);
    (*vec)[M+2*i] = *pos;
    (*vec)[M+2*i + 1] = *pos + 1;
    pool.erase(pos);

  }

}







void matrix_exp(vector<vector<double> > *mat, double m) {
	 
  // calculate a matrix exponential for FD conditionals
	 
  int i, j, n;
  double m1, m2;

  n = (*mat).size();
	 
  Matrix P(n, n), P1(n, n), B(n, n);

  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
			
      if(i == j) {
	B(i + 1, j + 1) = 1.0;
      }
      else {
	B(i + 1, j + 1) = 0.0;
      }

      P(i + 1, j + 1) = (*mat)[i][j];
			
    }
  }

  P1 = P;
  m1 = m;
  m2 = 1;
  for(i = 1; i < 100; i++) {
    B += (m1 / m2) * P1;
    P1 *= P;
    m1 *= m;
    m2 *= (i + 1);
  }

  B *= exp(-m);

  //cout << setw(10) << setprecision(5) << B << endl;
	
  for(i = 0; i < n; i++) 
    for(j = 0; j < n; j++) 
      (*mat)[i][j] = B(i + 1, j + 1);

}


int countLines (istream& in)
{
  return count(istreambuf_iterator<char>(in), istreambuf_iterator<char>(), '\n');
}


int rmultinom_unif(int n){

	/* generate a single draw from a Multinomial(probs) distribution */
	int i, ans = 0;
	double unif;
	Uniform runif;
	unif = runif.Next();
	i = 0;
	while(unif > 0 && ans < n) {
		ans = (i + 1);
		unif -= 1.0 / ((double) n);
		i += 1;
	}
	return ans;
}


int rmultinom(vector<double> *probs){

  //generate a single draw from a Multinomial(probs) distribution 
  int i = 0;
  //cout << "rmultinom" << endl;
  double unif;
  Uniform runif;
  unif = runif.Next();
  //cout << unif << endl;
  while(unif > 0 && i < (*probs).size()) 
    unif -= (*probs)[i++];
  //cout << i << endl;
  return i;
}


multimap<double, vector<int>, greater<double> >::iterator rmultinom(multimap<double, vector<int>, greater<double> > *m1, double total){

  // generate a single draw from a Multinomial(probs) distribution 
  int i = 0;
  double unif;
  Uniform runif;
  unif = runif.Next();
	
  multimap<double, vector<int>, greater<double> >::iterator pos;
  pos = (*m1).begin();

  while(unif > 0 && pos != (*m1).end()) {
    //cout << (pos->first / total) << endl;
    unif -= (pos->first / total);
    pos++;
  }
  pos--;
  return pos;
}


double genbet(float aa, float bb) {

  Gamma G1(aa), G2(bb);
  double g1, g2, beta;

  g1 = G1.Next();
  g2 = G2.Next();

  beta = g1 / (g1 + g2);

  return beta;
}


double rgamma(double shape, double rate) {

  // uses scaling property for Gamma distribution

  Gamma G(shape);
  double gam;

  gam = G.Next() / rate;

  return gam;

}

double genchi(int df) {
  
  double x;
  ChiSq chs(df);
  x = chs.Next();
  return x;
}


  double dgamma(double x, double a, double b, int flag) {
    
    /* calculates the gamma density */
    /* if flag = 1 it returns the log-density */
    /* if flag = 0 it returns the density */
    /* takes shape(=a) and rate(=b) */
    
    /* Density is given by f(x)= (b^a / Gamma(a)) x^(a-1) e^(-bx) */
    
    double d;
    
    d = a * log(b) + (b - 1) * x - b * x - lgamma(a);
    
    if(flag == 0) 
      return exp(d);
    else
      return d;
    
  }




void sequential_mean(double *mu, double *x, double m) {

	/* calculate a mean sequentially */

	*mu = (*x + (m - 1) * *mu) / m;

}

void sequential_variance(double *var, double *mu, double *x, double m) {

	/* calculate the variance sequentially */

	if(m == 1) *var = *x;
	if(m == 2) *var = ((*var - *x) / 2.0) * (*var - *x);
	if(m > 2) *var = (m - 2.0) * *var / (m - 1.0) + m * (*x - *mu) * (*x - *mu) / ((m - 1.0) * (m - 1.0)); 
}

void sequential_sd(double *sd, double *mu, double *x, double m) {

	/* calculate the variance sequentially */

	if(m == 1) *sd = *x;
	if(m == 2) *sd = sqrt(((*sd * *sd - *x) / 2.0) * (*sd * *sd - *x));
	if(m > 2) *sd = sqrt((m - 2.0) * *sd * *sd / (m - 1.0) + m * (*x - *mu) * (*x - *mu) / ((m - 1.0) * (m - 1.0))); 
}


void sequential_mean_vec(double *mu, double *x, double m, int vec_length) {

	int i;

	for(i = 0; i < vec_length; i++) sequential_mean(mu + i, x + i, m);
}

void sequential_variance_vec(double *var, double *mu, double *x, double m, int vec_length) {

	int i;

	for(i = 0; i < vec_length; i++) sequential_variance(var + i, mu + i, x + i, m); 
}

void sequential_sd_vec(double *sd, double *mu, double *x, double m, int vec_length) {

	int i;

	for(i = 0; i < vec_length; i++) sequential_sd(sd + i, mu + i, x + i, m); 
}
