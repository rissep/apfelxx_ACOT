//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/grid.h"
#include "apfel/dglap.h"
#include "apfel/observable.h"
#include "apfel/disbasis.h"

#include <functional>
#include <vector>

using std::function;
using std::vector;

namespace apfel
{
  /**
   * @brief Structure that contains all the precomputed quantities
   * needed to compute the DIS structure functions, i.e. the
   * perturbative coefficients of the coefficient functions for
   * F<SUB>2</SUB>, F<SUB>L</SUB>, and xF<SUB>3</SUB>.
   */
  struct StructureFunctionObjects
  {
    vector<int>             skip;
    map<int,ConvolutionMap> ConvBasis;
    map<int,Set<Operator>>  C0;
    map<int,Set<Operator>>  C1;
    map<int,Set<Operator>>  C2;
  };

  /**
   * @name DIS structure function object initializers
   * Collection of functions that initialise StructureFunctionObjects
   * structure for the different kinds of structure functions
   * available.
   */
  ///@{
  /**
   * @brief The InitializeF2NCObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for NC F2 in the ZM scheme
   * and store them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF2NCObjectsZM(Grid           const& g,
												   vector<double> const& Thresholds,
												   double         const& IntEps = 1e-5);

  /**
   * @brief The InitializeFLNCObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for NC FL in the ZM scheme
   * and store them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeFLNCObjectsZM(Grid           const& g,
												   vector<double> const& Thresholds,
												   double         const& IntEps = 1e-5);

  /**
   * @brief The InitializeF3NCObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for NC xF3 in the ZM scheme
   * and store them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF3NCObjectsZM(Grid           const& g,
												   vector<double> const& Thresholds,
												   double         const& IntEps = 1e-5);

  /**
   * @brief The InitializeF2CCPlusObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for combination ( F2(nu) +
   * F2(nubar) ) / 2 in the ZM scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF2CCPlusObjectsZM(Grid           const& g,
												       vector<double> const& Thresholds,
												       double         const& IntEps = 1e-5);

  /**
   * @brief The InitializeF2CCMinusObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for combination ( F2(nu) -
   * F2(nubar) ) / 2 in the ZM scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF2CCMinusObjectsZM(Grid           const& g,
													vector<double> const& Thresholds,
													double         const& IntEps = 1e-5);

  /**
   * @brief The InitializeFLCCPlusObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for combination ( FL(nu) +
   * FL(nubar) ) / 2 in the ZM scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeFLCCPlusObjectsZM(Grid           const& g,
												       vector<double> const& Thresholds,
												       double         const& IntEps = 1e-5);

  /**
   * @brief The InitializeFLCCMinusObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for combination ( FL(nu) -
   * FL(nubar) ) / 2 in the ZM scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeFLCCMinusObjectsZM(Grid           const& g,
													vector<double> const& Thresholds,
													double         const& IntEps = 1e-5);

  /**
   * @brief The InitializeF3CCPlusObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for combination ( F3(nu) +
   * F3(nubar) ) / 2 in the ZM scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF3CCPlusObjectsZM(Grid           const& g,
												       vector<double> const& Thresholds,
												       double         const& IntEps = 1e-5);

  /**
   * @brief The InitializeF3CCMinusObjectsZM precomputes the perturbative
   * coefficients of coefficient functions for combination ( F3(nu) -
   * F3(nubar) ) / 2 in the ZM scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF3CCMinusObjectsZM(Grid           const& g,
													vector<double> const& Thresholds,
													double         const& IntEps = 1e-5);

  /**
   * @brief The InitializeF2NCObjectsMassive precomputes the
   * perturbative coefficients of coefficient functions for
   * combination NC F2 in the massive scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nxi: the number of nodes of the grid in &xi; = Q<SUP>2</SUP>/M<SUP>2</SUP> (default: 150)
   * @param ximin: the lower bound of the grid in &xi; (default: 0.001)
   * @param ximax: the upper bound of the grid in &xi; (default: 100000)
   * @param intdeg: the interpolation degree on the grid in &xi; (default: 3)
   * @param lambda: the value of the parameter in the function ln(ln(&xi;/&Lambda;<SUP>2</SUP>)) used for the tabulation (default: 0.0005)
   * @return A StructureFunctionObjects-valued function
   */
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF2NCObjectsMassive(Grid           const& g,
													vector<double> const& Masses,
													double         const& IntEps = 1e-5,
													int            const& nxi    = 150,
													double         const& ximin  = 0.001,
													double         const& ximax  = 100000,
													int            const& intdeg = 3,
													double         const& lambda = 0.0005);

  /**
   * @brief The InitializeFLNCObjectsMassive precomputes the
   * perturbative coefficients of coefficient functions for
   * combination NC FL in the massive scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nxi: the number of nodes of the grid in &xi; = Q<SUP>2</SUP>/M<SUP>2</SUP> (default: 150)
   * @param ximin: the lower bound of the grid in &xi; (default: 0.001)
   * @param ximax: the upper bound of the grid in &xi; (default: 100000)
   * @param intdeg: the interpolation degree on the grid in &xi; (default: 3)
   * @param lambda: the value of the parameter in the function ln(ln(&xi;/&Lambda;<SUP>2</SUP>)) used for the tabulation (default: 0.0005)
   * @return A StructureFunctionObjects-valued function
   */
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeFLNCObjectsMassive(Grid           const& g,
													vector<double> const& Masses,
													double         const& IntEps = 1e-5,
													int            const& nxi    = 150,
													double         const& ximin  = 0.001,
													double         const& ximax  = 100000,
													int            const& intdeg = 3,
													double         const& lambda = 0.0005);

  /**
   * @brief The InitializeF2NCObjectsMassiveZero precomputes the
   * perturbative coefficients of coefficient functions for
   * combination NC F2 in the massless limit of the massive scheme and
   * store them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nxi: the number of nodes of the grid in &xi; = Q<SUP>2</SUP>/M<SUP>2</SUP> (default: 150)
   * @param ximin: the lower bound of the grid in &xi; (default: 0.001)
   * @param ximax: the upper bound of the grid in &xi; (default: 100000)
   * @param intdeg: the interpolation degree on the grid in &xi; (default: 3)
   * @param lambda: the value of the parameter in the function ln(ln(&xi;/&Lambda;<SUP>2</SUP>)) used for the tabulation (default: 0.0005)
   * @return A StructureFunctionObjects-valued function
   */
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeF2NCObjectsMassiveZero(Grid           const& g,
													    vector<double> const& Masses,
													    double         const& IntEps = 1e-5,
													    int            const& nxi    = 150,
													    double         const& ximin  = 0.001,
													    double         const& ximax  = 100000,
													    int            const& intdeg = 3,
													    double         const& lambda = 0.0005);

  /**
   * @brief The InitializeFLNCObjectsMassiveZero precomputes the
   * perturbative coefficients of coefficient functions for
   * combination NC FL in the massless limit of the massive scheme and
   * store them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nxi: the number of nodes of the grid in &xi; = Q<SUP>2</SUP>/M<SUP>2</SUP> (default: 150)
   * @param ximin: the lower bound of the grid in &xi; (default: 0.001)
   * @param ximax: the upper bound of the grid in &xi; (default: 100000)
   * @param intdeg: the interpolation degree on the grid in &xi; (default: 3)
   * @param lambda: the value of the parameter in the function ln(ln(&xi;/&Lambda;<SUP>2</SUP>)) used for the tabulation (default: 0.0005)
   * @return A StructureFunctionObjects-valued function
   */
  function<StructureFunctionObjects(double const&, vector<double> const&)> InitializeFLNCObjectsMassiveZero(Grid           const& g,
													    vector<double> const& Masses,
													    double         const& IntEps = 1e-5,
													    int            const& nxi    = 150,
													    double         const& ximin  = 0.001,
													    double         const& ximax  = 100000,
													    int            const& intdeg = 3,
													    double         const& lambda = 0.0005);
  ///@}

  /**
   * @name Structure function builders
   * Collection of functions that build a map of Observable objects
   * corresponding to the different component of the structure
   * functions.
   */
  ///@{
  /**
   * @brief The StructureFunctionBuildNC class constructs a map of
   * "Observable" objects.
   * @param FObj: the StructureFunctionObjects-valued for the structure function objects
   * @param InDistFunc: the distribution to be convoluted with as a map<int,double>-valued function of x and Q
   * @param PerturbativeOrder: the perturbative order
   * @param Alphas: the strong coupling function
   * @param Couplings: the vector-valued function of (non-QCD) couplings
   * @return a 
   */
  map<int,Observable> BuildStructureFunctions(function<StructureFunctionObjects(double const&, vector<double> const&)> const& FObj,
					      function<map<int,double>(double const&, double const&)>                  const& InDistFunc,
					      int                                                                      const& PerturbativeOrder,
					      function<double(double const&)>                                          const& Alphas,
					      function<vector<double>(double const&)>                                  const& Couplings);

  /**
   * @brief The StructureFunctionBuildNC class constructs a map of
   * "Observable" objects.
   * @param FObj: the StructureFunctionObjects-valued for the structure function objects
   * @param InDistFunc: the distribution to be convoluted with as a double-valued function of i, x, and Q
   * @param PerturbativeOrder: the perturbative order
   * @param Alphas: the strong coupling function
   * @param Couplings: the vector-valued function of (non-QCD) couplings
   * @return a 
   */
  //_____________________________________________________________________________
  map<int,Observable> BuildStructureFunctions(function<StructureFunctionObjects(double const&, vector<double> const&)> const& FObj,
					      function<double(int const&, double const&, double const&)>               const& InDistFunc,
					      int                                                                      const& PerturbativeOrder,
					      function<double(double const&)>                                          const& Alphas,
					      function<vector<double>(double const&)>                                  const& Couplings);
  ///@}
}
