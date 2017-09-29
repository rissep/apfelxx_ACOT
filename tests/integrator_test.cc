//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/integrator.h>

#include <iostream>
#include <cmath>

using namespace apfel;
using namespace std;

int main()
{ 
  const Integrator f{[&] (double const& x)->double{ return log(x); }};
  // Integrate using dgauss with a given accuracy.
  const double res1 = f.integrate(0, 2, 1e-5);
  // Integrate using dgauss with a given number of nodes.
  const double res2 = f.integrate(0, 2, 4);

  cout << res1 << "  " << res2 << "  " << res1 / res2 << endl;

  // Integrate introducing fixed point in the integration interval.
  const double res3 = f.integrate(0, 2, {-1, 0.1, 0.5, 0.7, 1.1, 3.4}, 1e-5);

  cout << res1 << "  " << res3 << "  " << res1 / res3 << endl;

  // Now revert integration and try integrate with and without fixed
  // points.
  const double res4 = f.integrate(2, 0, 1e-5);
  const double res5 = f.integrate(2, 0, {-1, 0.1, 0.5, 0.7, 1.1, 3.4}, 1e-5);

  cout << res4 << "  " << res5 << "  " << res4 / res5 << endl;

  return 0;
}
