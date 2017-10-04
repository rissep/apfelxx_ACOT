//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <functional>
#include <vector>
#include <iostream>

using std::function;
using std::vector;

namespace apfel
{
  /**
   * @brief The Integrator class which uses Adaptative Gaussian Quadrature.
   *
   * This class takes as input the integrand function and provides the
   * integrate method which performs the integration.
   */
  class Integrator
  {
  public:
    /**
     * @brief The default constructor
     */
    Integrator();

    /**
     * @brief The default constructor
     * @param func the function to be integrated.
     */
    Integrator(function<double(double const&)> const& func);

    /**
     * @brief This takes a function of two variables and integrates
     * over the first keeping the second fixed to 'arg2'.
     */
    Integrator(function<double(double const&, double const&)> const& func2, double const& arg2);

    /**
     * @brief This takes a function of three variables and integrates
     * over the first keeping the second and the thisrd fixed to
     * 'arg2' and 'arg3'.
     */
    Integrator(function<double(double const&, double const&, double const&)> const& func3, double const& arg2, double const& arg3);

    /**
     * @brief Integrates the integrand passed during initialization
     * between xmin and xmax with tolerance eps.
     *
     * @param xmin the lower bound integration value.
     * @param xmax the upper bound integration value.
     * @param eps the required relative error.
     * @return the integral value.
     */
    double integrate(double const& xmin, double const& xmax, double const& eps) const;

    /**
     * @brief Integrates the integrand passed during initialization
     * between xmin and xmax with tolerance eps.
     *
     * @param xmin the lower bound integration value.
     * @param xmax the upper bound integration value.
     * @param FixPts vector of fixed points not to integrate over.
     * @param eps the required relative error.
     * @return the integral value.
     */
    double integrate(double const& xmin, double const& xmax, vector<double> const& FixPts, double const& eps) const;

    /**
     * @brief Integrates the integrand passed during initialization
     * between xmin and xmax with m points.
     *
     * @param xmin the lower bound integration value.
     * @param xmax the upper bound integration value.
     * @param m the number of points of the Gauss quadrature.
     * @return the integral value.
     */
    double integrate(double const& xmin, double const& xmax, int const& m) const;

  protected:
    /**
     * @brief Protected virtual integrand function.
     *
     * @param x the integration variable.
     * @return the integrand evaluated at x.
     */
    virtual double integrand(double const& x) const { return _func(x); };

  private:
    function<double(double const&)> _func;
  };
}
