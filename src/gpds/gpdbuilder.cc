//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gpdbuilder.h"
#include "apfel/dglap.h"
#include "apfel/operator.h"
#include "apfel/operatorerbl.h"
#include "apfel/set.h"
#include "apfel/timer.h"
#include "apfel/constants.h"
#include "apfel/tools.h"
#include "apfel/messages.h"
#include "apfel/gpdevolutionkernels.h"
#include "apfel/evolutionbasisqcd.h"
#include "apfel/matchingbasisqcd.h"

namespace apfel
{
  //_____________________________________________________________________________
  std::map<int, DglapObjects> InitializeGpdObjects(Grid                const& g,
                                                   std::vector<double> const& Thresholds,
                                                   double              const& xi,
                                                   bool                const& OpEvol,
                                                   double              const& IntEps)
  {
    report("Initializing DglapObjects for GPD unpolarised evolution... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the thresholds
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // Allocate needed operators (matching conditions and splitting
    // functions). By now the code is fast enough to precompute
    // everything at all available perturbative orders and the current
    // perturbative order is accounted for only when the actual
    // splitting functions and matching conditions (lambda) functions
    // are defined.
    // ===============================================================
    // LO Matching conditions.
    std::map<int, Operator> MatchLO;
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};
    MatchLO.insert({MatchingBasisQCD::M0, Id});
    MatchLO.insert({MatchingBasisQCD::M1, Zero});
    MatchLO.insert({MatchingBasisQCD::M2, Zero});
    MatchLO.insert({MatchingBasisQCD::M3, Zero});
    MatchLO.insert({MatchingBasisQCD::M4, Zero});
    MatchLO.insert({MatchingBasisQCD::M5, Zero});
    MatchLO.insert({MatchingBasisQCD::M6, Zero});
    MatchLO.insert({MatchingBasisQCD::M7, Zero});

    // ===============================================================
    // LO splitting function operators.
    std::map<int, std::map<int, Operator>> OpMapLO;
    const Operator O0ns = Operator{g, Pgpd0nsDGLAP{xi}, IntEps};// + OperatorERBL{g, Pgpd0nsERBL{xi}, IntEps};
    const Operator O0qg = Operator{g, Pgpd0qgDGLAP{xi}, IntEps};// + OperatorERBL{g, Pgpd0qgERBL{xi}, IntEps};
    const Operator O0gq = Operator{g, Pgpd0gqDGLAP{xi}, IntEps};// + OperatorERBL{g, Pgpd0gqERBL{xi}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O0gg = Operator{g, Pgpd0ggDGLAP{nf, xi}, IntEps} + OperatorERBL{g, Pgpd0ggERBL{nf, xi}, IntEps};
        const Operator O0qgnf = nf * O0qg;
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O0ns});
        OM.insert({EvolutionBasisQCD::PNSM, O0ns});
        OM.insert({EvolutionBasisQCD::PNSV, O0ns});
        OM.insert({EvolutionBasisQCD::PQQ,  O0ns});
        OM.insert({EvolutionBasisQCD::PQG,  O0qgnf});
        OM.insert({EvolutionBasisQCD::PGQ,  O0gq});
        OM.insert({EvolutionBasisQCD::PGG,  O0gg});
        OpMapLO.insert({nf, OM});
      }

    // Define object of the structure containing the DglapObjects.
    std::map<int, DglapObjects> DglapObj;

    // Allocate convolution maps for evolution and matching, and set
    // of operators.
    for (int nf = nfi; nf <= nff; nf++)
      {
        DglapObjects obj;
        obj.Threshold = Thresholds[nf-1];
        if (OpEvol)
          {
            obj.SplittingFunctions.insert({0, Set<Operator>{EvolutionOperatorBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.MatchingConditions.insert({0, Set<Operator>{MatchingOperatorBasisQCD{nf},  MatchLO}});
          }
        else
          {
            obj.SplittingFunctions.insert({0, Set<Operator>{EvolutionBasisQCD{nf}, OpMapLO.at(nf)}});
            obj.MatchingConditions.insert({0, Set<Operator>{MatchingBasisQCD{nf},  MatchLO}});
          }
        DglapObj.insert({nf,obj});
      }
    t.stop();

    return DglapObj;
  }
}
