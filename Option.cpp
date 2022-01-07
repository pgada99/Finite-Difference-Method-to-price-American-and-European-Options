// Pietro Gadaleta / Francesco Stampini
// Finite Differences Option Pricing
// C++ Project

/*! \file option.cpp contains definitions of methods for the two classes defined in option.h:
     - mesh:    used to define a time mesh grid and a mesh grid for the spot
                price (using derived class "stock_mesh")
     - option:  used to define EU/AM calls and puts. the option contains all the
                methods necessary to compute the option's price and greeks.*/

#include <iostream>
#include "Option.h"
/*!namespace where are defined all the classes used*/
namespace m2qf{




//! mesh constructor - used for time mesh.
mesh::mesh(double p_0, double p_fin, unsigned long N): N_(N), dt_((p_fin-p_0)/N), p0_(p_0)
{}




//! stock mesh constructor.
stock_mesh::stock_mesh(double p_0, double p_fin, unsigned long N, double p): mesh(p_0, p_fin, N)
{
    // p = value of stock at time 0. we want it to be
    // included in the mesh.
    //grid_.reserve(N);
    // creating grid s.t. p is in the grid.

    unsigned long i=1;
    dt_ = (p-p_0)/i;
    while((p_0+dt_*N)>p_fin){
        i++;
        dt_ = (p-p_0)/i;
    }
    
    idx_= i;  // storing S0's index in grid_.
    //for(unsigned long i=0; i<N+1; ++i)
    //    grid_.push_back(p_0+i*dt_);
}




//! The tridiagonal matrix algorithm (TDMA), also known als Thomas algorithm, is a simplified form of Gaussian elimination that can be used to solve tridiagonal system of equations.
void option::TDMA_Algorithm(const std::vector<double>& a, const std::vector<double>& b,
                            const std::vector<double>& c, const std::vector<double>& y,
                            std::vector<double>& x, std::vector<double>& c_tmp,
                            std::vector<double>& y_tmp)
{
   unsigned long m = y.size(); // m = M-1
   c_tmp[0] = c[0]/b[0];
   y_tmp[0] = y[0]/b[0];

   double n;
   for (unsigned long i=1; i<m; ++i)
   {
       n = 1.0/(b[i]-a[i-1]*c_tmp[i-1]); //is a[i-1] since we have erased the first term of the vector a
       c_tmp[i] = c[i]*n;
       y_tmp[i] = (y[i]-a[i-1]*y_tmp[i-1])*n;
   }
   // reverse sweep to update solution vector:
   for (unsigned long i = m-1; i-- > 0; )
   {
       if(i==m-1)
           x[i] = y_tmp[i]-c_tmp[i]*y_tmp[i+1];
       else
           x[i] = y_tmp[i]-c_tmp[i]*x[i+1];
   }
}




//! Payoff function
double get_payoff (bool isCall, double S, double K)
{
    return (isCall)*std::max(S-K,0.0) + (!isCall)*std::max(K-S,0.0);
}





//! Function to price the option
double option::price_solve( double r, double sigma, const stock_mesh & mesh_spot, const mesh & mesh_time, bool isCall, bool isAmerican)
{
    std::vector<double> a, b, c; // coefficients of tridiagonal matrix
    std::vector<double> y;       // right-hand column vector
    std::vector<double> result;  // vector giving option values for each spot value in the mesh
                                 // (eventually, it will return those values at t=0)
    
    // tridiagonal matrix is an (M-1)x(M-1) matrix.
    // so, at the end:
    // a.size() = M-2
    // b.size() = M-1
    // c.size() = M-2
    
    // just for now we will be making all diagonals of size M-1.
    // we will shrink a and c in a later step.
    a.reserve(mesh_spot.get_size()-2);
    b.reserve(mesh_spot.get_size()-2);
    c.reserve(mesh_spot.get_size()-2);

    // y.size() = M-1 (same size as tridiagonal matrix)
    y.reserve(mesh_spot.get_size()-2);

    // result.size() should be M+1 (same size as mesh).
    // But we are removing boundary values, since we know that, at every timepoint, we have
    // f(i,0) = K and f(i,M) = 0 (for a put option, for all i = 0,...,N).
    result.resize(mesh_spot.get_size()-2);
    
    double dt = mesh_time.get_dt();
    double ds = mesh_spot.get_dt();
    long N = mesh_time.get_size()-1;
    long M = mesh_spot.get_size()-1;
    
    if (isAmerican) // resizing ex_b_tmp
        ex_b_tmp.resize(mesh_time.get_size());
    
    // creating the tridiagonal matrix with the coefficients (see formulas in the report):
    for(unsigned long j=1; j < mesh_spot.get_size()-1; ++j) // j = 1:M
    {
        a.push_back(0.5*dt/ds*(r-pow(sigma,2)/2)-pow(sigma/ds, 2)*0.5*dt);
        b.push_back(1 + pow(sigma/ds, 2)*dt + r*dt);
        c.push_back(-pow(sigma/ds, 2)*0.5*dt -0.5*(r-pow(sigma,2)/2)*dt/ds);
        
        // starting at time t=T (maturity): value of option is equal to payoff (K-S)+
        y.push_back(get_payoff(isCall, exp(mesh_spot[j]), strike_)); // y now contains f(i+1,j) for i=N+1 (t=T) and j=1,...,M-1
    }
    
    // storing boundary values of a and c that will not be in the tridiagonal matrix.
    // (i.e. initial value of a(j=1) and final value of c(j=M+1))
    // we will need them (one or the other, depending on option type) to put in the right-hand
    // column vector.
    double a0=a[0];
    double cM1 = c[mesh_spot.get_size()-3]; //c(M-1)

    if (isCall)
    {
        //  y[end] ~ f(i+1,M) - c(j=M+1)*(S_max - K*exp(T-t))
        y[mesh_spot.get_size()-3] -= cM1*(exp(mesh_spot.get_max())-strike_);
    }
    else // isPut
    {
        // removing initial value of a(j=1). before that, we will store it.
        // in fact, y[0] = f(i+1,1) - a(j=1)*K  for all i=0...N-1, since f(i,0) = K for all i.
        y[0] -= a0 * strike_;
    }

    // removing values
    a.erase(a.begin());
    c.pop_back();

    // creating two vectors here that will be used for the TDMA algorithm function.
    std::vector<double> ctmp(y.size()), ytmp(y.size());
    
    for(unsigned long n=1; n<mesh_time.get_size(); ++n)
    {
        
        //initialising c_tmp, y_tmp
        for (unsigned long i = 0; i < ytmp.size(); ++i)
        {
            ctmp[i]=0.0;
            ytmp[i]=0.0;
        }
        
        TDMA_Algorithm(a,b,c,y,result,ctmp,ytmp);
        ctmp.clear();
        ytmp.clear();
       
        
         // adding continuity condition if option is american
         if (isAmerican)
         {
             bool ok = true;
             // default value of ex_b
             ex_b_tmp[N-n] = std::exp(mesh_spot_[mesh_spot_.get_size()-1]);
             if(isCall)
             {
                 for (long int i = 0; i < result.size(); ++i)
                 {
                  // if condition above holds we are in the boundary between continuation region and exercise region: we store index in the ex_b_tmp
                     if(result[i]<get_payoff(isCall, std::exp(mesh_spot[i+1]), strike_) && ok==true)
                     {
                         ok = false;
                         
                         ex_b_tmp[N-n] = std::exp(mesh_spot_[i]);//mesh spot is a step bigger than result
                     }
                     result[i] = std::max(result[i], get_payoff(isCall, std::exp(mesh_spot[i+1]), strike_));
                 }
             }
             if(!isCall)
             {
                 for (long int i = result.size()-1; i >=0 ; --i)
                 {
                     // if condition above holds we are in the boundary between continuation region and exercise region: we store index in the ex_b_tmp
                     if(result[i]<get_payoff(isCall, std::exp(mesh_spot[i+1]), strike_) && ok==true)
                     {
                         ok = false;
                         
                         ex_b_tmp[N-n] = std::exp(mesh_spot_[i]);//mesh spot is a step bigger than result
                     }
                     result[i] = std::max(result[i], get_payoff(isCall, std::exp(mesh_spot[i+1]), strike_));
                 }
             }
         }
         // end of american condition


        // if   y      = [ f(i+1,1)-a(j=1)*K, f(i+1,2), ..., f(i+1,M-1)]
        // then result = [ f(i,1), ..., f(i,M-1) ]

        // in order to update it for the next iteration we set:
        y.clear();
        y = result;
        
        // ...and...
        if (isCall)
            y[mesh_spot.get_size()-3] -= cM1*(exp(mesh_spot.get_max())-strike_*exp(-r*mesh_time.get_dt()*n));
        else // isPut
            y[0] -= a0*strike_;
    };
    
    // if the option is American, we are adding extreme values to exercise boundary vector
    if(isAmerican)
    {
        //ex_b_tmp[mesh_time_.get_size()-1] = std::exp(mesh_spot_.get_max());
        ex_b_tmp[0] = ex_b_tmp[1];
    }
    // the value of the call/put option will therefore be positioned in the vector
    // "result" at the stored index. Remember to subtract one since
    // result starts from j=1 and not j=0!
    int pos = mesh_spot.get_index()-1;
    double result__ = result[pos];
    
    price_vector_tmp = result; // storing results for each spot's price in the spot mesh grid.
    return result__;
}

    
//! function used to compute the greek of the option + plots
void option::greek_solve(bool isCall, bool isAmerican)
{
    
    // the following vectors are needed to plot graphs of greeks vs. underlying spot price.
    // when we reach the actual spot price, we store result for the option we are considering.
    stock_vect_.reserve(interval_+1);
    price_vect_.reserve(interval_+1);
    delta_vect_.reserve(interval_+1);
    if (extra_plots_)
    {
        theta_vect_.reserve(interval_+1);
        rho_vect_.reserve(interval_+1);
        vega_vect_.reserve(interval_+1);
    }
    
    int s_idx = mesh_spot_.get_index();
    int cidx;// current index)
    int dM= std::floor((mesh_spot_.get_size()-1)/(2*interval_)); // step of the grid
    double h = 0.001; //step for volatility, time, interest rate
    double vegaprice=0, rhoprice=0, thetaprice=0, ctrprice=0;
                      // prices we will need to calculate for estimating greeks theta, rho, vega.
    
    // this for-loop is in order to build all vectors for plots.
    // when we are reaching the actual spot price, we store result for the option we are considering.
    for (int i = 0; i <= interval_; ++i)
    {
            //cidx = s_idx-std::floor(3*interval_/4)+i;     // defining current index
            cidx = s_idx-std::floor(dM*interval_/2)+i*dM;
            
            stock_vect_.push_back(exp(mesh_spot_[cidx]));
            price_vect_.push_back(price_vector_[cidx]);
            // we use the center bump formula for derivative
            delta_vect_.push_back((price_vector_[cidx+1]-price_vector_[cidx-1])/(exp(mesh_spot_[cidx+1])-exp(mesh_spot_[cidx-1])));
            //gamma price is calculated using the second derivative formula approximation
            //f(x)'' = (f(x+h)-2*f(x)+f(x-h ))/h^2
            gamma_vect_.push_back((price_vector_[cidx+1]+price_vector_[cidx-1]-2*price_vector_[cidx])/((exp(mesh_spot_[cidx+1])-exp(mesh_spot_[cidx]))*(exp(mesh_spot_[cidx])-exp(mesh_spot_[cidx-1])))); // assuming the two ds (mesh spots) are similar.
            
            // defining stock_mesh
            if(extra_plots_ || cidx == s_idx) // if we want plots for vega, theta, rho OR
                                              // we need to compute greeks for our option
            {
                stock_mesh mesh_spot_center(mesh_spot_[cidx]-l(sigma_,mesh_time_.get_max(),r_), mesh_spot_[cidx]+l(sigma_,mesh_time_.get_max(),r_), mesh_spot_.get_size()-1, mesh_spot_[cidx]);
                // computing prices for greeks (vega, rho, theta) by adding increment h to derivation variable
                ctrprice = price_solve(r_, sigma_, mesh_spot_center, mesh_time_, isCall, isAmerican);
                vegaprice = price_solve(r_, sigma_+h, mesh_spot_center, mesh_time_,isCall, isAmerican);
                rhoprice = price_solve(r_+h, sigma_, mesh_spot_center,mesh_time_,isCall, isAmerican);
                mesh mesh_temp(0.0,mesh_time_.get_max()-h,mesh_time_.get_size()-1);
                thetaprice = price_solve(r_, sigma_,mesh_spot_, mesh_temp, isCall, isAmerican);
            }
            if(extra_plots_)
            {
                // inserting estimation of greeks into vectors
                vega_vect_.push_back((vegaprice-ctrprice)/h);
                rho_vect_.push_back((rhoprice-ctrprice)/(h));
                theta_vect_.push_back((thetaprice-ctrprice)/h);
            }
            if (cidx == s_idx) // if current idx is equal to spot price's index: storing values of greeks
            {
                delta_  = delta_vect_[i];//(delta_vect_[i]+delta_vect_[i-1])/2;
                gamma_  = gamma_vect_[i];//(gamma_vect_[i]+gamma_vect_[i-1])/2;
                vega_   = (vegaprice-ctrprice)/h;
                rho_    = (rhoprice-ctrprice)/h;
                theta_  = (thetaprice-ctrprice)/h;
            }
    }
}

    
    
    
//!function used to return the scalar value of the greek
double option::get_greek(const std::string & greek) const
{
    if(greek=="delta")
        return delta_;
    else if(greek=="rho")
        return rho_;
    else if(greek=="vega")
        return vega_;
    else if(greek=="theta")
        return theta_;
    else if(greek=="gamma")
        return gamma_;
    else
        throw "EXCEPTION";
}

    
    
    
//!function used to return the vector of the price or greek in the interval
std::vector<double> option::get_vect_plot(const std::string & name)const{
    if(name=="delta")
        return delta_vect_;
    else if(name=="price")
        return price_vect_;
    else if(name=="gamma")
        return gamma_vect_;
    else if(name=="vega")
        return vega_vect_;
    else if(name=="rho")
        return rho_vect_;
    else if(name=="theta")
        return theta_vect_;
    else if(name=="stock")
        return stock_vect_;
    else
        throw "EXCEPTION";
}
    
}// end ofm2qf namespace

