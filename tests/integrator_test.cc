//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <iostream>
#include <cmath>

int main()
{
  std::cout << "- 1D integration" << std::endl;

  // Integrand function.
  const apfel::Integrator fGK{[&] (double const& x) -> double{ return log(x); }};
  const apfel::Integrator fGL{[&] (double const& x) -> double{ return log(x); }, apfel::Integrator::IntegrationMethod::GAUSS_LEGENDRE};

  // Print true value.
  std::cout << "True value: " << 2 * ( log(2) - 1 ) << std::endl;

  // Integrate using Gauss-Legendre with a given accuracy.
  const double res1 = fGL.integrate(0, 2, 1e-5);

  // Integrate using Gauss-Legendre with a different accruracy.
  const double res2 = fGL.integrate(0, 2, 1e-3);

  // Print results.
  std::cout << "Gauss-Legendre: " << res1 << "  " << res2 << "  " << res1 / res2 << std::endl;

  // Integrate using Gauss-Kronrod with a given accuracy.
  const double res6 = fGK.integrate(0, 2, 1e-5);

  // Integrate using Gauss-Kronrod with a different accruracy.
  const double res7 = fGK.integrate(0, 2, 1e-3);

  // Print results.
  std::cout << "Gauss-Kronrod:  " << res6 << "  " << res7 << "  " << res6 / res7 << std::endl;

  // Now revert integration and try integrate with and without fixed
  // points.
  const double res4 = fGK.integrate(0, 2, 1e-5);
  const double res5 = fGK.integrate(0, 2, {-1, 0.1, 0.5, 0.7, 1.1, 3.4}, 1e-5);

  // Print results.
  std::cout << "W/o and w/ fixed points: " << res4 << "  " << res5 << "  " << res4 / res5 << std::endl;

  // Performance test
  int k = 10000;
  std::cout << "Integrating " << k << " times with Gauss-Legendre... ";
  apfel::Timer t;
  for (int i = 0; i < k; i++)
    fGL.integrate(0, 2, 1e-9);
  t.stop();

  std::cout << "Integrating " << k << " times with Gauss-Kronrod... ";
  t.start();
  for (int i = 0; i < k; i++)
    fGK.integrate(0, 2, 1e-9);
  t.stop();

  // 2D integration
  std::cout << "\n- 2D integration" << std::endl;

  // Integrand function.
  const apfel::Integrator2D f2d{[&] (double const& x, double const& y) -> double{ return log(x) * log(y); }};

  // Print true value.
  std::cout << "True value: " << 4 * pow(log(2) - 1, 2) << std::endl;

  // Integrate using Gauss-Legendre with a given accuracy.
  const double res2d1 = f2d.integrate(0, 2, 0, 2, 1e-7);

  // Integrate using Gauss-Legendre with a different accruracy.
  const double res2d2 = f2d.integrate(0, 2, 0, 2, 1e-3);

  // Print results.
  std::cout << "Gauss-Legendre: " << res2d1 << "  " << res2d2 << "  " << res2d1 / res2d2 << std::endl;

  // Performance test
  k = 100;
  std::cout << "Integrating " << k << " times with Gauss-Legendre... ";
  t.start();
  for (int i = 0; i < k; i++)
    f2d.integrate(0, 2, 0, 2, 1e-7);
  t.stop();

  return 0;
}
