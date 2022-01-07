// Finite Differences Option Pricing
// C++ Project

/*!\file toexcel.cpp contains the implementation of the function "get_results" which we pass to Excel VBA. */

#include <stdio.h>
#include "toexcel.h"
#include "Option.h"


double get_results(double strike, double S_0, double r, double sigma, double t_fin, double Nd, double Md, double calld, double americand, double intervald, double extra_plotsd, double * results, double * price, double * greek_vect)//, double * gamma_vect)//, double * rho_vect, double * vega_vect)
{
    // call and american are double, either equal to 0.0 or 1.0
    bool isCall = calld;
    bool isAmerican = americand;
    bool extraplots = extra_plotsd; // if we want to add plots for rho, theta, vega too.
    
    long N = Nd;
    long M = Md;
    
    //interval of the graph
    int interval = intervald;
    // Nd, Md are integers but of type double, we want to make them long integers
        
    // data is pointing an array where i am going to put price + greeks
    m2qf::option op(strike,S_0,r,sigma,t_fin,
                    N,M,isCall,isAmerican, interval,extraplots);
    
    results[0] = op.get_result();
    results[1] = op.get_greek("delta");
    results[2] = op.get_greek("gamma");
    results[3] = op.get_greek("rho");
    results[4] = op.get_greek("theta");
    results[5] = op.get_greek("vega");
    
    std::vector<double> price_vector= op.get_vect_plot("price");
    // greek is initialised with the array of stock_prices for plots
    std::vector<double> greek_vector= op.get_vect_plot("stock");
    std::vector<double> delta_vector= op.get_vect_plot("delta");
    std::vector<double> gamma_vector= op.get_vect_plot("gamma");
    // in case of extraplots
    std::vector<double> rho_vector= op.get_vect_plot("rho");
    std::vector<double> vega_vector= op.get_vect_plot("vega");
    std::vector<double> theta_vector= op.get_vect_plot("theta");
    
    // appending greeks arrays to greek_vector
    // reserving space
    if(extraplots)
        greek_vector.reserve((interval+1)*6);
    else
        greek_vector.reserve((interval+1)*3);
    
    greek_vector.insert(greek_vector.end(), delta_vector.begin(), delta_vector.end());
    greek_vector.insert(greek_vector.end(), gamma_vector.begin(), gamma_vector.end());
    if (extraplots)
    {
        // adding rho, vega, theta vectors for plots
        std::vector<double> rho_vector= op.get_vect_plot("rho");
        std::vector<double> vega_vector= op.get_vect_plot("vega");
        std::vector<double> theta_vector= op.get_vect_plot("theta");
        
        greek_vector.insert(greek_vector.end(), rho_vector.begin(), rho_vector.end());
        greek_vector.insert(greek_vector.end(), vega_vector.begin(), vega_vector.end());
        greek_vector.insert(greek_vector.end(), theta_vector.begin(), theta_vector.end());
    }
    
    // copying the result in two arrays that we pass by reference to VBA
    for(int i=0; i<=interval; i++)
        price[i]= price_vector[i];
    
    // vector "price" will also contain the values to plot the exercise boundary
    if(isAmerican)
    {
        for (int i=0; i<=N; i++)
            price[interval+1+i] = op.get_ex_b(i);
        
    }
    
    for(int i=0; i<greek_vector.size();i++)
                greek_vect[i]=greek_vector[i];
    
    return op.get_result();
}
