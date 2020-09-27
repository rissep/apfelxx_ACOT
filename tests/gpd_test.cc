//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <functional>

int main()
{
  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100,1e-5,3}, apfel::SubGrid{60,1e-1,3}, apfel::SubGrid{50,6e-1,3}, apfel::SubGrid{50,8e-1,3}}};

  // Initial scale
  const double mu0 = sqrt(2);

  // Vector of thresholds
  const std::vector<double> Thresholds = {0, 0, 0, sqrt(2), 4.5, 175};

  // Perturbative order
  const int PerturbativeOrder = 0;

  // Skewness
  const double xi = 0.01;

  // Running coupling
  const double AlphaQCDRef = 0.35;
  const double MuAlphaQCDRef = sqrt(2);
  apfel::AlphaQCD a{AlphaQCDRef, MuAlphaQCDRef, Thresholds, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // Initialize GPD evolution objects
  const auto GpdObj   = InitializeGpdObjects(g, Thresholds, xi);
  const auto GpdObjOp = InitializeGpdObjects(g, Thresholds, xi, true);

  // Construct the DGLAP objects
  const auto EvolvedGPDs = BuildDglap(GpdObj,   apfel::LHToyPDFs, mu0, PerturbativeOrder, as);
  const auto EvolvedOps  = BuildDglap(GpdObjOp,                   mu0, PerturbativeOrder, as);

  // Tabulate GPDs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedGPDs{*EvolvedGPDs, 50, 1, 1000, 3};

  // Tabulate Operators
  const apfel::TabulateObject<apfel::Set<apfel::Operator>> TabulatedOps{*EvolvedOps, 50, 1, 1000, 3};

  // Final scale
  const double mu = 10;

  // Compute results
  std::cout << std::scientific << "Direct evolution (4th order Runge-Kutta) from Q0 = " << mu0 << " GeV to Q = " << mu << " GeV... ";

  // Evolve GPDs to the final Scale
  apfel::Timer t;
  const std::map<int, apfel::Distribution> gpds = apfel::QCDEvToPhys(EvolvedGPDs->Evaluate(mu).GetObjects());
  t.stop();

  std::cout << "Interpolation of the tabulated evolution operators... ";
  t.start();
  apfel::Set<apfel::Operator> tops = TabulatedOps.Evaluate(mu);
  t.stop();

  // Set appropriate convolution basis for the evolution operators and
  // convolute them with initial-scale distributions.
  tops.SetMap(apfel::EvolveDistributionsBasisQCD{});
  const std::map<int, apfel::Distribution> opgpds = apfel::QCDEvToPhys((tops * apfel::Set<apfel::Distribution> {apfel::EvolveDistributionsBasisQCD{}, DistributionMap(g, apfel::LHToyPDFs, mu0)}).GetObjects());

  // Print results
  const std::vector<double> xlha = {1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

  std::cout << "\nAlphaQCD(Q) = " << Alphas.Evaluate(mu) << std::endl;
  std::cout << "\n   x    "
            << "   u-ubar   "
            << "   d-dbar   "
            << " 2(ubr+dbr) "
            << "   c+cbar   "
            << "    gluon   "
            << std::endl;

  std::cout << "Direct Evolution:" << std::endl;
  for (auto const& x : xlha)
    {
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << "  " << (gpds.at(2) - gpds.at(-2)).Evaluate(x)
                << "  " << (gpds.at(1) - gpds.at(-1)).Evaluate(x)
                << "  " << 2 * (gpds.at(-2) + gpds.at(-1)).Evaluate(x)
                << "  " << (gpds.at(4) + gpds.at(-4)).Evaluate(x)
                << "  " << gpds.at(0).Evaluate(x)
                << std::endl;
    }
  std::cout << "\n";

  std::cout << "Evolution through the evolution operator:" << std::endl;
  for (auto const& x : xlha)
    {
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << "  " << (opgpds.at(2) - opgpds.at(-2)).Evaluate(x)
                << "  " << (opgpds.at(1) - opgpds.at(-1)).Evaluate(x)
                << "  " << 2 * (opgpds.at(-2) + opgpds.at(-1)).Evaluate(x)
                << "  " << (opgpds.at(4) + opgpds.at(-4)).Evaluate(x)
                << "  " << opgpds.at(0).Evaluate(x)
                << std::endl;
    }
  std::cout << "\n";

  return 0;
}
