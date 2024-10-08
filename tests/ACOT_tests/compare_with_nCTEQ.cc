#include <fstream>
#include <iomanip>

#include "apfel/apfelxx.h"
#include "apfel/structurefunctionbuilderACOT.h"
#include "LHAPDF/LHAPDF.h"

const std::vector<double> x{
  1.00000000e-06, 1.14975700e-06, 1.32194115e-06, 1.51991108e-06,
  1.74752840e-06, 2.00923300e-06, 2.31012970e-06, 2.65608778e-06,
  3.05385551e-06, 3.51119173e-06, 4.03701726e-06, 4.64158883e-06,
  5.33669923e-06, 6.13590727e-06, 7.05480231e-06, 8.11130831e-06,
  9.32603347e-06, 1.07226722e-05, 1.23284674e-05, 1.41747416e-05,
  1.62975083e-05, 1.87381742e-05, 2.15443469e-05, 2.47707636e-05,
  2.84803587e-05, 3.27454916e-05, 3.76493581e-05, 4.32876128e-05,
  4.97702356e-05, 5.72236766e-05, 6.57933225e-05, 7.56463328e-05,
  8.69749003e-05, 1.00000000e-04, 1.14975700e-04, 1.32194115e-04,
  1.51991108e-04, 1.74752840e-04, 2.00923300e-04, 2.31012970e-04,
  2.65608778e-04, 3.05385551e-04, 3.51119173e-04, 4.03701726e-04,
  4.64158883e-04, 5.33669923e-04, 6.13590727e-04, 7.05480231e-04,
  8.11130831e-04, 9.32603347e-04, 1.07226722e-03, 1.23284674e-03,
  1.41747416e-03, 1.62975083e-03, 1.87381742e-03, 2.15443469e-03,
  2.47707636e-03, 2.84803587e-03, 3.27454916e-03, 3.76493581e-03,
  4.32876128e-03, 4.97702356e-03, 5.72236766e-03, 6.57933225e-03,
  7.56463328e-03, 8.69749003e-03, 1.00000000e-02, 1.14975700e-02,
  1.32194115e-02, 1.51991108e-02, 1.74752840e-02, 2.00923300e-02,
  2.31012970e-02, 2.65608778e-02, 3.05385551e-02, 3.51119173e-02,
  4.03701726e-02, 4.64158883e-02, 5.33669923e-02, 6.13590727e-02,
  7.05480231e-02, 8.11130831e-02, 9.32603347e-02, 1.07226722e-01,
  1.23284674e-01, 1.41747416e-01, 1.62975083e-01, 1.87381742e-01,
  2.15443469e-01, 2.47707636e-01, 2.84803587e-01, 3.27454916e-01,
  3.76493581e-01, 4.32876128e-01, 4.97702356e-01, 5.72236766e-01,
  6.57933225e-01, 7.56463328e-01, 8.69749003e-01
};
const std::vector<double> Q{1.4,5,10,100};

int main(){

  LHAPDF::PDF* pdf = LHAPDF::mkPDF("CT18NNLO",0);
  const auto PDFrotated = [&] (double const& x, double const& Q) -> std::map<int,double>{return apfel::PhysToQCDEv(pdf->xfxQ(x,Q));};
  const std::vector<double> Thresholds = {0,0,0,pdf->quarkMass(4),pdf->quarkMass(5),pdf->quarkMass(6)};
  const auto alphas = [&] (double const& Q) -> double{return pdf->alphasQ(Q);};

  double IntEps = 1e-5;
  int nQ = 60;
  double Qmin = 1.29;
  double Qmax = 101;
  int intdeg = 3;
  int pto = 1;
  apfel::Grid g{{{25,1e-6,4},{20,1e-2,4},{10,1e-1,4},{5,5e-1,4}}};

  const auto fEW = [=] (double const& Q) -> std::vector<double> {return apfel::ElectroWeakCharges(Q,false);};

  // const auto F2objects = apfel::InitializeF2NCObjectsMassive(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg);
  const auto F2objects = apfel::InitializeF2NCObjectsACOT(g,Thresholds,IntEps,nQ,Qmin,Qmax,intdeg);
  const std::map<int,apfel::Observable<>> F2 = apfel::BuildStructureFunctions(F2objects,PDFrotated,pto,alphas,fEW);
  const apfel::TabulateObject<apfel::Distribution> F2total {[&] (double const& Q) -> apfel::Distribution{return F2.at(0).Evaluate(Q);},nQ,Qmin,Qmax,intdeg,Thresholds};

  std::string path = "/home/users/p_riss01/codes/apfelxx_ACOT/tests/ACOT_tests/APFELxx_results";
  std::string filename = path + "/NC_F2_ACOT.csv";

  std::ofstream file;
  file.open(filename);
  file<<std::fixed<<std::setprecision(10);
  file<<"x,Q2,NCF2\n";
  for(double xi: x){
    for(double Qi: Q){
      file<<xi<<","<<Qi*Qi<<",";
      file<<F2total.EvaluatexQ(xi,Qi);
      file<<"\n"; 
    }
  }
  file.close();

  return 0;
}