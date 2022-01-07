// Pietro Gadaleta / Francesco Stampini
// Finite Differences Option Pricing
// C++ Project

/*! \file option.h contains definitions of two classes:
     - mesh:    used to define a time mesh grid and a mesh grid for the spot
                price (using derived class "stock_mesh")
     - option:  used to define EU/AM calls and puts. the option contains all the
                methods necessary to compute the option's price and greeks.*/
#ifndef OPTION_H
#define OPTION_H

/* The classes below are exported */
#pragma GCC visibility push(default)

#include <stdio.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
/*!namespace where are defined all the classes used*/
namespace m2qf{

//! Class defining a regular mesh given starting and final point and length of the grid.
class mesh{
    
    protected:
        //! number of points = N+1
        unsigned long       N_;
        //! delta between mesh values
        double              dt_;
        //! initial point
        double              p0_;
    
    public:
        mesh(){}
        //!Constructor
        mesh(double p_0, double p_fin, unsigned long N);
        //! return the point in that index of the grid
        double operator[](size_t i) const {return p0_ + i*dt_;}
        //! return the size of the grid
        unsigned long get_size() const { return  N_+1;}
        //! return the difference between two point of the grid
        double get_dt() const{return dt_;}
        //! return the final point of the grid
        double get_max() const {return p0_ + N_*dt_;}
};

//! Class used to create a mesh for the stock's price.
/*! This class is made in order to build a mesh for the stock prices.
    Compared to its mother class "mesh", it needs an extra element, i.e. the index of the place in which S_0 is stored in the grid. */
class stock_mesh: public mesh

{
    unsigned long idx_; // index of the mesh indicating S0's position.
    public:
        //!Constructor
        stock_mesh(double p_0, double p_fin, unsigned long N, double p);
        double get_index() const{return  idx_;}
        
};

//! Option class used to price an European or American option.
/*! This class is made in order to compute the price and the greeks of the option. The option can be a simple European option (Call or Put), or an American option (Call or Put). */
class option
{
    protected:
        // data
        double strike_; // strike price
        double r_;      // risk-free rate
        double sigma_;  // volatility
        double S_0_;    // underlying price at t = 0 (today)
        bool isCall_;   // = 1 if option is Call;
                        // = 0 if option is Put
        bool isAmerican_; // = 1 if option is American;
                        //   = 0 if option is European
        
        mesh mesh_time_;        // time mesh class
        stock_mesh mesh_spot_;  // underlying spot mesh class
        
        // results:
        double price_;  // price of the option
        double delta_;  // greeks: delta (dP/dS)
        double gamma_;  // greeks: gamma (d2P/dS2)
        double rho_;    // greeks: rho   (dP/dr)
        double vega_;   // greeks: vega  (dP/dsigma)
        double theta_;  // greeks: theta (-dP/dt)
    
        //! utilities:
        //! boolean: if true, we also compute prices for vega, theta and rho plots. otherwise, only delta and gamma are computed for plotting. this is to cut down computational speed
        bool extra_plots_;
        // vector results
        //! vector of option's prices for each underlying stock's price in the mesh grid.
        std::vector<double> price_vector_;
        //! vector of option's  greek and price for each underlying stock's price in the mesh grid to return to VBA for the plot.
        std::vector<double> delta_vect_, gamma_vect_, price_vect_, rho_vect_,
                            vega_vect_, theta_vect_, stock_vect_;
        //! this vector changes everytime price_solve is called. it changes from price_vector_ as soon as we call greek_solve, since we are changing the parameters of price_solve.
        std::vector<double> price_vector_tmp;
        //! vector containing indices of the stock_mesh for each time point to indicate the exercise boundary for an American option.
        std::vector<double> ex_b;
        //!this vector changes everytime price_solve is called. it changes from ex_b as soon as we call greek_solve, since we are changing the parameters of price_solve.
        std::vector<double> ex_b_tmp;
        //! number of points of the grid to return to VBA for the plots
        int interval_;
    
    protected:
        //######FUNCTION
        
        //!The tridiagonal matrix algorithm (TDMA), also known als Thomas algorithm, is a simplified form of Gaussian elimination that can be used to solve tridiagonal system of equations.
        /*!This function solves the problem Ax = y, where A is a tridiagonal matrix.
         Inputs:
             a,b,c = vectors forming tri-diagonal matrix A. b is the main diagonal, a is the lower second-diagonal and
                     c is the upper one.
             y     = right-hand side vector. same size as matrix A.
             x     = solution vector. it will be overwritten.
             c_tmp
             y_tmp = vectors used in the algorithm. since we are going to run this function N times, we
                     initialize them outside of the function to save excess allocation/deallocation.
                     they have size (M-1) and are initialized by zero.
         */
        void TDMA_Algorithm(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c,
                            const std::vector<double>& y, std::vector<double>& x,
                            std::vector<double>& c_tmp, std::vector<double>& y_tmp);
    
        //! Function that computes the boundaries of the stock mesh grid, l returns the semi-length of the spot mesh interval
        double l(double sigma, double t_fin, double r)
        {
            return abs(r-pow(sigma,2)/2)*t_fin + 6*sigma*pow(t_fin,0.5);
        }
    
        //! Function that computes the price of the option
        double price_solve(double r, double sigma, const stock_mesh & mesh_spot, const mesh & mesh_time, bool isCall = false, bool isAmerican = false);
    
        //! Function that computes the payoff of the option
        void payoff(double S_T, double K, const std::string & option);
    
        //! Function that computes all the greeks of the option
        void greek_solve(bool isCall = false, bool isAmerican = false);
        
    public:
        //!Constructor of the class option
        option(double strike, double S_0, double r, double sigma, double t_fin, unsigned long  N, unsigned long M, bool isCall, bool isAmerican, int interval, bool extra_plots):
                        strike_(strike), r_(r), sigma_(sigma),
                        mesh_spot_(log(S_0)-l(sigma,t_fin,r), log(S_0)+l(sigma,t_fin,r), M, log(S_0)),
                        S_0_(S_0), interval_(interval),
                        mesh_time_(0.0, t_fin, N),
                        isCall_(isCall), isAmerican_(isAmerican), extra_plots_(extra_plots)
        {
            price_ = price_solve(r_, sigma_, mesh_spot_, mesh_time_, isCall_, isAmerican_);
            price_vector_ = price_vector_tmp; // values of option's price are stored in price_vector.
                                                  // we store it in a different vector since price_vector_tmp
                                                  // changes everytime we call price_solve in greek_solve.
            ex_b = ex_b_tmp;
            greek_solve(isCall_, isAmerican_);
        }
        
        //! function that returns the price
        double get_result() const {return price_;}
        
        //!function that returns greeks, please specify the name of the greek (case sensitive)
        double get_greek(const std::string & greek) const;
        
        //! function to get values for the exercise boundary
        double get_ex_b(const long i) const {return ex_b[i];}
        
        //!function that returns the vector to plot, so price and the vector of each greek. Returns a std vector.
        std::vector<double> get_vect_plot(const std::string & name)const;
    
}; // end of class option

} // end of namespace m2qf

#pragma GCC visibility pop
#endif //OPTION_H
