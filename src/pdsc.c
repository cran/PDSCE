#include <math.h>
#include <stdlib.h>


double **make_mat(int nrow, int ncol);
void delete_mat( double **mat);
double *make_vec(int len);
void delete_vec(double *vec);
double soft(double x, double t);                  
				  
void pdsc(double * Sin, double *Cov0, double * Inv0, int * pin, 
            double * lamin, 
            double * bin,
            double * tolin, 
            int * maxitin, 
            double * tolout, 
            int * maxitout, 
            int * totalout, 
            double * Sout,
            double * Oout)
{
  int p=*pin;
  double b=*bin;
  double **S=make_mat(p,p);
  double **Sigma=make_mat(p,p);
  double **Omega=make_mat(p,p); 
  double **lam=make_mat(p,p);
  
  //  make p rows and only one column
  double *beta=make_vec(p); 
  //used for the lasso regression
  double *coltotal=make_vec(p);  
  /* used for the termination criteria for the lasso regression */
  
  int i,j,jj,k, numit_in, numit_out;
  double diff_in, diff_out, stotal, tmp;
  double ct;
  //  Initialize iterates and read in the sample cov/cor
  stotal=0;
  for(i=0; i < p; i++)
  {
    ct=0;
    for ( j=0; j < p; j++)
    {
      S[i][j]=Sin[j*p+i];
      lam[i][j]=lamin[j*p+i];
      Sigma[i][j]=Cov0[j*p+i];
      Omega[i][j]=Inv0[j*p+i];  
      if( i < j)
        stotal += fabs(S[i][j]);
      if( i != j )
        ct=ct+fabs(S[i][j]);
    }
    coltotal[i]=ct;
  }

  numit_out=0;
  diff_out=*tolout*stotal+1;
  while( (diff_out > (*tolout * stotal) ) && (numit_out < *maxitout) )
  {
    numit_out+=1;
    diff_out = 0;
    for(j=0; j < p; j++)
    {
      // update jth diagonal entry of Sigma
      tmp = S[j][j] + b*Omega[j][j]; 
      diff_out += fabs(tmp - Sigma[j][j]);
      Sigma[j][j] = tmp;
       
      numit_in=0;
     
      //initialize beta    
      for(k=0; k < p; k++)
      {      
        beta[k] = Sigma[k][j]; 
      } 

      diff_in=*tolin*coltotal[j]+1;      
      while( (diff_in > (*tolin*coltotal[j])) && (numit_in < *maxitin) )
      {
        numit_in+=1;
        diff_in = 0;
        for(jj=0; jj < p; jj++)
        {        
          if(jj != j)
          {
            tmp = 0;
            for(k=0; k<p; k++)
            {
              if(k != jj && k != j)
              {
               tmp=tmp+ Omega[k][jj] * beta[k];
              }
            }
            tmp = S[jj][j] - b*tmp/Sigma[j][j];
            tmp = soft(tmp, lam[j][jj]);
            tmp = tmp/( 1 + b*Omega[jj][jj]/Sigma[j][j]);
            diff_in += fabs(tmp - beta[jj]);
            beta[jj] = tmp;
          }
        }
      } // end inner while

      /* the lasso regression for the jth row/column is complete 
         with the result stored in beta */      

      // move beta into the jth row/column of Sigma (no diagonal entry)
      for(jj=0; jj <p; jj++)
      {
        if(jj != j)
        {
          tmp = beta[jj];
          diff_out += fabs(tmp - Sigma[jj][j]);
          Sigma[jj][j] = tmp;
          Sigma[j][jj] = tmp;
        }
      }

      // update the jth row/column of Omega (no diagonal entry)
      for(jj=0; jj <p; jj++)
      {
        if(jj != j)
        {
          tmp=0;
          for(k=0; k < p; k++)
          {
            if(k != j)
              tmp = tmp + Omega[jj][k] * beta[k];
          }
          tmp=-tmp/Sigma[j][j];
          Omega[jj][j] = tmp;
          Omega[j][jj] = tmp;
        }
      }

      // update Omega jj
      tmp=0;
      for(k=0; k < p; k++)
      {
        if(k != j)
          tmp = tmp + beta[k] * Omega[k][j];
      }
      tmp = 1 - tmp;
      tmp = tmp/Sigma[j][j];
      Omega[j][j] = tmp;
    }// end for loop over j   
  } // end outer while
  totalout[0] = numit_out;
  //prepare output into Bout
  for(j=0; j <p; j++)
  {
    for (i=0; i < p; i++)
    {    
      Sout[j*p+i] = Sigma[i][j];
      Oout[j*p+i] = Omega[i][j];
    }
  }
  delete_mat(S);
  delete_mat(Sigma);
  delete_mat(Omega);
  delete_mat(lam);
  delete_vec(beta);
  delete_vec(coltotal);
}


				  
void bchol(double * xin, int * nin, int * pin, int * kin, double * bcov)
{
  int p=*pin;
  int n=*nin;
  int k=*kin;
  int i,ii,iii,j,kk,q;
  double xtx,xty,rss,tmpsum,sum,tmp,tmp2,d,a,b;
  double **L=make_mat(p,p); 
  /* L is the covariance cholesky factor, 
    it's diagonal is D */
  double **x=make_mat(n,p);
  double **r=make_mat(n,p); //matrix of residuals
  
  // Read in input 
  for(j=0; j < p; j++)
  {
    for ( i=0; i < n; i++)
    {
      x[i][j]=xin[j*n+i];
      r[i][j]=xin[j*n+i];
    }
  }
  
  // get first residual variance
  rss=0;
  for(ii=0;ii < n; ii++)
  {
    rss+=x[ii][0]*x[ii][0];
  }  
  L[0][0] = rss/n;
  
  
  // go thru rows of cholesky factor (starting with the second row)  
  for( i=1; i <p; i++)
  {
    //get coefficients in the ith row of L
    for(j=1; (j<=k) && ((i-j) >= 0); j++)
    {
      xtx=0;
      xty=0;
      for(ii=0; ii < n; ii++)
      {
        xtx += r[ii][i-j] *r[ii][i-j];
        xty += r[ii][i-j] *x[ii][i];
      }
      L[i][i-j] = xty/xtx; 
    }
    //get residuals
    rss=0;
    for(ii=0; ii < n; ii++)
    {
      sum=0;
      for(iii=1; (iii <= k) && ((i-iii) >= 0); iii++)
      {
        sum += r[ii][i-iii] * L[i][i-iii];
      }
      r[ii][i] = x[ii][i]-sum;
      rss += r[ii][i]*r[ii][i];	   
    }
    L[i][i]=rss/n;
  }
  
  for(i=1; i <=p; i++)
  {
    for(j = i; (i-j) <= k && j >= 1; j--)
    {
      tmpsum=0;
      for(kk=0; (kk <= k) && (kk <= (k-i+j)) &&((j-kk) >= 1); kk++)
      {
        //for cov i,j
	q = j-kk;
	d = L[q-1][q-1];
	a = L[i-1][q-1];
	b = L[j-1][q-1];
	if ( i == q)
	{
	  a=1;
	}
	if( j == q)
	{
	  b=1;
	}
        tmpsum += a * d * b;
      }
      bcov[(j-1)*p+(i-1)] = tmpsum;
      bcov[(i-1)*p+(j-1)] = tmpsum;
    }
  }
  delete_mat(L);
  delete_mat(x);
  delete_mat(r);
}

double soft(double x, double t)
{
  double tmp, ax, newx;
  newx=0;
  ax = fabs(x);
  tmp = ax - t;
  if( tmp > 0)
  {
    if(x > 0)
    {
      newx = tmp;
    }
    else
    {
      newx = -tmp;
    }
  }
  return newx;
}


double **make_mat(int nrow, int ncol)
{
  double ** mat;
  int k;
  
  mat = (double **) malloc(nrow*sizeof(double*));
  mat[0]=(double*) malloc(nrow*ncol*sizeof(double));
  
  for(k=1; k < nrow; k++)
    mat[k] = mat[k-1] + ncol;
  return mat;
}
void delete_mat( double **mat)
{
  free(mat[0]);
  free(mat);  
}


double *make_vec(int len)
{
  double * vec;
  vec = (double *) malloc(len*sizeof(double));
  return vec;
}

void delete_vec(double *vec)
{
  free(vec);  
}


