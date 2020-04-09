//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/lagrangeinterpolator.h"
#include "apfel/constants.h"

#include <algorithm>

namespace apfel
{
  //_________________________________________________________________________________
  LagrangeInterpolator::LagrangeInterpolator(Grid const& gr):
    Interpolator{gr}
  {
  }

  //_________________________________________________________________________________
  LagrangeInterpolator::LagrangeInterpolator(Grid                             const& gr,
                                             std::vector<std::vector<double>> const& distsubgrid,
                                             std::vector<double>              const& distjointgrid):
    Interpolator{gr, distsubgrid, distjointgrid}
  {
  }

  //_________________________________________________________________________________
  double LagrangeInterpolator::Interpolant(int const& beta, double const& lnx, SubGrid const& sg) const
  {
    // Get the logarithmic grid.
    const std::vector<double>& lxsg = sg.GetLogGrid();

    // Return immediately 1 if "x" coincides with "xg[beta]".
    if (std::abs(lnx - lxsg[beta]) < eps12)
      return 1;

    // Define the lower bound of the interpolation range.
    const int id    = sg.InterDegree();
    const int bound = std::max(beta - id, 0);

    // Return 0 if "x" is outside the range in which the interpolant
    // is different from zero.  Ideally this functions should never be
    // called if "beta" and "x" are such that "Interpolant" is
    // identically zero. Use "SumBounds" to know where "beta" should
    // run over given "x".
    if (lnx < lxsg[bound] || lnx >= lxsg[beta+1])
      return 0;

    // Find the the neighbors of "x" on the grid.
    int j;
    for (j = 0; j <= beta - bound; j++)
      if (lnx >= lxsg[beta-j])
        break;

    // Compute the interpolant.
    double w_int = 1;
    for (int delta = beta - j; delta <= beta - j + id; delta++)
      if (delta != beta)
        w_int *= ( lnx - lxsg[delta] ) / ( lxsg[beta] - lxsg[delta] );

    return w_int;
  }

  //_________________________________________________________________________________
  double LagrangeInterpolator::DerInterpolant(int const& beta, double const& lnx, SubGrid const& sg) const
  {
    // Get the logarithmic grid.
    const std::vector<double>& lxsg = sg.GetLogGrid();

    // Define the lower bound of the interpolation range.
    const int id    = sg.InterDegree();
    const int bound = std::max(beta - id, 0);

    // Return 0 if "x" is outside the range in which the interpolant
    // is different from zero.  Ideally this functions should never be
    // called if "beta" and "x" are such that "Interpolant" is
    // identically zero. Use "SumBounds" to know where "beta" should
    // run over given "x".
    if (lnx < lxsg[bound] || lnx >= lxsg[beta+1])
      return 0;

    // Find the the neighbors of "x" on the grid.
    int j;
    for (j = 0; j <= beta - bound; j++)
      if (lnx >= lxsg[beta-j])
        break;

    // Compute the interpolant.
    double dw_int = 0;
    for (int gamma = beta - j; gamma <= beta - j + id; gamma++)
      {
        double w = 1;
        for (int delta = beta - j; delta <= beta - j + id; delta++)
          if (delta != beta && delta != gamma)
            w *= ( lnx - lxsg[delta] ) / ( lxsg[beta] - lxsg[delta] );
        if (gamma != beta)
          {
            w /= lxsg[beta] - lxsg[gamma];
            dw_int += w;
          }
      }
    return dw_int;
  }

  //_________________________________________________________________________________
  std::array<int, 2> LagrangeInterpolator::SumBounds(double const& x, SubGrid const& sg) const
  {
    const std::vector<double>& xsg = sg.GetGrid();

    std::array<int,2> bounds = {{0, 0}};
    if (x < xsg[0] - eps12 || x > xsg[sg.nx()] + eps12)
      return bounds;

    const int low = std::lower_bound(xsg.begin(), xsg.end() - sg.InterDegree() - 1, x) - xsg.begin();
    bounds[0] = bounds[1] = low;

    if (std::abs(x - xsg[low]) <= eps12)
      bounds[1] += 1;
    else
      {
        bounds[0] -= 1;
        bounds[1] += sg.InterDegree();
      }

    return bounds;
  }
}
