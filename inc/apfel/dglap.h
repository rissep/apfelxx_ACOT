//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/matchedevolution.h"
#include "apfel/set.h"
#include "apfel/distribution.h"
#include "apfel/operator.h"

#include <functional>

using std::function;

namespace apfel
{
  /**
   * @brief The Dglap class.
   *
   * A specialization class of the MatchedEvolution class
   * for the computation of the DGLAP evolution.
   */
  class Dglap: public MatchedEvolution<Set<Distribution>>
  {
  public:

    Dglap() =  delete;

    /**
     * @brief Dglap default constructor.
     */
    Dglap(function<Set<Operator>(int,double)>      const& SplittingFunctions,
	  function<Set<Operator>(bool,int,double)> const& MatchingConditions,
	  Set<Distribution>                        const& ObjRef,
	  double                                   const& MuRef,
	  vector<double>                           const& Masses,
	  vector<double>                           const& Thresholds,
	  int                                      const& nsteps = 10);

    /**
     * @brief Dglap
     */
    Dglap(function<Set<Operator>(int,double)>      const& SplittingFunctions,
	  function<Set<Operator>(bool,int,double)> const& MatchingConditions,
	  Set<Distribution>                        const& ObjRef,
	  double                                   const& MuDistRef,
	  vector<double>                           const& Masses,
	  int                                      const& nsteps = 10);

    /**
     * @brief
     * @param
     * @return
     */
    Set<Distribution> MatchObject(bool const& Up, int const& nf, Set<Distribution> const& sd) const;

    /**
     * @brief
     * @param
     * @return
     */
    Set<Distribution> Derivative(int const& nf, double const& mu, Set<Distribution> const& f) const;

    /**
     * @brief Function that returns the number of steps.
     */
    int const& GetNumberOfSteps() const { return _nsteps; }

  private:
    function<Set<Operator>(int,double)>      _SplittingFunctions;
    function<Set<Operator>(bool,int,double)> _MatchingConditions;
    int                               const  _nsteps;
  };

}
