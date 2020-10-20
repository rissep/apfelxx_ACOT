//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/grid.h"
#include "apfel/expression.h"
#include "apfel/distribution.h"

namespace apfel
{
  /**
   * @brief The Operator class defines the basic object "Operator"
   * which is essentially the convolution on the grid bewteen an
   * Expression object (e.g. a splitting function) and the interpolant
   * functions.
   */
  class Operator
  {
  public:
    Operator() = delete;

    /**
     * @brief The Operator constructor.
     * @param gr: the Grid object
     * @param expr: the expression to be transformed
     * @param eps: relative accuracy of the numerical integrations (default: 10<SUP>-5</SUP>)
     * @param erbl: whether the convolution integral is ERBL like (default: false)
     */
    Operator(Grid const& gr, Expression const& expr, double const& eps = 1e-5, bool const& erbl = false);

    /**
     * @brief The Operator constructor.
     * @param gr: the Grid object
     * @param expr: the expression to be transformed
     * @param erbl: whether the convolution integral is ERBL like
     * @param eps: relative accuracy of the numerical integrations (default: 10<SUP>-5</SUP>)
     */
    Operator(Grid const& gr, Expression const& expr, bool const& erbl, double const& eps = 1e-5);

    /**
     * @brief Function that returns the ERBL-like flag
     */
    bool IsERBL() const { return _erbl; }

    /**
     * @name Binary operators
     */
    ///@{
    Distribution operator *= (Distribution const& d) const;            //!< this *= Distribution
    Operator&    operator  = (Operator const& o);                      //!< this  = Operator
    Operator&    operator *= (Operator const& o);                      //!< this *= Operator
    Operator&    operator *= (double const& s);                        //!< this *= Scalar
    Operator&    operator *= (std::function<double(double const&)> f); //!< This *= Function
    Operator&    operator /= (double const& s);                        //!< this /= Scalar
    Operator&    operator += (Operator const& o);                      //!< this += Operator
    Operator&    operator -= (Operator const& o);                      //!< this -= Operator
    ///@}

    /**
     * @brief Function that returns the Grid object of the operator.
     */
    Grid const& GetGrid() const { return _grid; }

    /**
     * @brief Function that returns the operator.
     */
    std::vector<std::vector<double>> GetOperator() const { return _Operator; }

  protected:
    Grid                        const& _grid;         //!< Grid on which to compute the operator
    bool                        const  _erbl;         //!< If the operator is ERBL-like
    std::vector<std::vector<double>>   _Operator;     //!< Operator values.

    friend std::ostream& operator << (std::ostream& os, Operator const& sg);
  };

  /**
   * @name Ternary operators
   */
  ///@{
  Distribution operator * (Operator lhs, Distribution const& rhs);                //!< Operator*Distribution
  Operator     operator * (Operator lhs, Operator const& rhs);                    //!< Operator*Operator
  Operator     operator * (double const& s, Operator rhs);                        //!< Scalar*Operator
  Operator     operator * (Operator lhs, double const& s);                        //!< Operator*Scalar
  Operator     operator * (std::function<double(double const&)> f, Operator rhs); //!< function*Operator
  Operator     operator * (Operator lhs, std::function<double(double const&)> f); //!< Operator*function
  Operator     operator / (Operator lhs, double const& s);                        //!< Operator/Scalar
  Operator     operator + (Operator lhs, Operator const& rhs);                    //!< Operator+Operator
  Operator     operator - (Operator lhs, Operator const& rhs);                    //!< Operator-Operator
  ///@}

  /**
   * @brief Method which prints Operator with cout <<. This only
   * prints the Operator on the first subgrid and is supposed to be
   * used for debugging purposes.
   */
  std::ostream& operator << (std::ostream& os, Operator const& op);
}
