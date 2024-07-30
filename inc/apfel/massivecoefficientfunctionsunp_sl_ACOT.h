#pragma once

#include "apfel/expression.h"

namespace apfel
{

  /**
   * @defgroup ACOT Massiv coefficient functions
   * Collection of massive coefficient functions for the
   * (full) ACOT scheme up to O((&alpha;<SUB>s</SUB>).
   * @note Only neutral current is available.
   */

  ///@{
  class Cm20qNC_ACOT: public Expression{
    public:
      Cm20qNC_ACOT(double const& eta);
      double Local(double const& x) const;

    private:
      double _xi;
  };

  class Cm21gNC_ACOT: public Expression
  {
  public:
    Cm21gNC_ACOT(double const& eta);
    double Regular(double const& x) const;
  };

  class Cm21gNC_sub_ACOT: public Expression
  {
  public:
    Cm21gNC_sub_ACOT(double const& eta);
    double Regular(double const& x) const;
  private:
    double _xi;
  };

  class Cm21qNC_ACOT: public Expression
  {
  public:
    Cm21qNC_ACOT(double const& eta);
    double Local(double const& x) const;
    double Singular(double const& x) const;
  private:
    double g(double const& x) const;
    double _sud;
    double _V2;
    double _S2;
    double _xi;
  };

  class Cm21qNC_sub_ACOT: public Expression
  {
  public:
    Cm21qNC_sub_ACOT(double const& eta);
    double Local(double const& x) const;
    double Singular(double const& x) const;
  private:
    double _xi;
  };

  class C31nsNC_ACOT : public Expression
  {
  public:
    C31nsNC_ACOT(double const& eta);
    double Singular(double const& x) const;
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
  };

  class CL1gNC_ACOT : public Expression
  {
  public:
    CL1gNC_ACOT(double const& eta);
    double Regular(double const& x)  const;
  };

  class CL1nsNC_ACOT : public Expression
  {
  public:
    CL1nsNC_ACOT(double const& eta);
    double Regular(double const& x)  const;
  };
  ///@}

  /**
   * @defgroup SACOT-chi Massiv coefficient functions
   * Collection of massive coefficient functions for the
   * SACOT-chi scheme up to O((&alpha;<SUB>s</SUB>).
   */
  ///@{
  class Cm20qNC_ACOT_chi: public Expression{
    public:
      Cm20qNC_ACOT_chi(double const& eta);
      double Local(double const& x) const;
  };

  class Cm21gNC_sub_ACOT_chi: public Expression
  {
  public:
    Cm21gNC_sub_ACOT_chi(double const& eta);
    double Regular(double const& x) const;
  };

  class Cm21qNC_ACOT_chi : public Expression
  {
  public:
    Cm21qNC_ACOT_chi(double const& eta);
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  class Cm21gCC_sub_ACOT: public Expression
  {
  public:
    Cm21gCC_sub_ACOT(double const& eta, double const& xi);
    double Regular(double const& x)  const;
  private:
    const double _xi;
  };

  class Cm21gCC_general_mass: public Expression
  {
  public:
    Cm21gCC_general_mass(double const& eta, double const& xi, double const& m1, double const& m2);
    double Regular(double const& x)  const;
  private:
    const double _xi;
    const double _m1;
    const double _m2;
  };

  //unify naming scheme?
  class CmL1gCC_general_mass: public Expression 
  {
  public:
    CmL1gCC_general_mass(double const& eta, double const& xi, double const& m1, double const& m2);
    double Regular(double const& x)  const;
  private:
    const double _xi;
    const double _m1;
    const double _m2;
  };

  class Cm31gCC_general_mass: public Expression 
  {
  public:
    Cm31gCC_general_mass(double const& eta, double const& xi, double const& m1, double const& m2);
    double Regular(double const& x)  const;
  private:
    const double _xi;
    const double _m1;
    const double _m2;
  };

  class C21CCg_one_mass: public Expression
  {
  public:
    C21CCg_one_mass(double const& eta, double const& xi);
    double Regular(double const& x)  const;
  private:
    const double _xi;
  };
  ///@}

  /**
   * @defgroup aSACOT-chi Massiv coefficient functions
   * Collection of massive coefficient functions for the
   * approximative SACOT-chi scheme only for 
   * O((&alpha;<SUB>s</SUB><SUP>2</SUP>).
   */
  ///@{
  class C22nsm_ACOT_NNLO_0 : public Expression
  {
  public:
    C22nsm_ACOT_NNLO_0(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };

  class C22nsm_ACOT_NNLO_nf : public Expression
  {
  public:
    C22nsm_ACOT_NNLO_nf(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };

  class C22nsp_ACOT_NNLO_0 : public Expression
  {
  public:
    C22nsp_ACOT_NNLO_0(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };

  class C22nsp_ACOT_NNLO_nf : public Expression
  {
  public:
    C22nsp_ACOT_NNLO_nf(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };

  class CL2nsm_ACOT_NNLO_0 : public Expression
  {
  public:
    CL2nsm_ACOT_NNLO_0(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
  };

  class CL2nsm_ACOT_NNLO_nf : public Expression
  {
  public:
    CL2nsm_ACOT_NNLO_nf(double const& eta,bool const& xdependent);
    double Regular(double const& x)  const;
  };

  class CL2g_ACOT_NNLO : public Expression
  {
  public:
    CL2g_ACOT_NNLO(double const& eta,bool const& xdependent);
    double Regular(double const& x)  const;
  };

  class CL2ps_ACOT_NNLO : public Expression
  {
  public:
    CL2ps_ACOT_NNLO(double const& eta,bool const& xdependent);
    double Regular(double const& x)  const;
  };

  class CL2nsp_ACOT_NNLO_0 : public Expression
  {
  public:
    CL2nsp_ACOT_NNLO_0(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
  };

  class CL2nsp_ACOT_NNLO_nf : public Expression
  {
  public:
    CL2nsp_ACOT_NNLO_nf(double const& eta,bool const& xdependent);
    double Regular(double const& x)  const;
  };

  class C32nsm_ACOT_NNLO_0 : public Expression
  {
  public:
    C32nsm_ACOT_NNLO_0(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };

  class C32nsm_ACOT_NNLO_nf : public Expression
  {
  public:
    C32nsm_ACOT_NNLO_nf(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };

  class C32nsp_ACOT_NNLO_0 : public Expression
  {
  public:
    C32nsp_ACOT_NNLO_0(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };

  class C32nsp_ACOT_NNLO_nf : public Expression
  {
  public:
    C32nsp_ACOT_NNLO_nf(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };
  ///@}




  

  

  

  //sim ACOT-chi NNLO trick
  //Error in name? should be C22g??
  class C21g_ACOT_NNLO: public Expression
  {
  public:
    C21g_ACOT_NNLO(double const& eta);
    double Regular(double const& x) const;
  };

  class C22nsp_ACOT_NNLO: public Expression
  {
  public:
    C22nsp_ACOT_NNLO(int const& nf, double const& eta, bool const& xdependent);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };

  class C22ps_ACOT_NNLO: public Expression
  {
  public:
    C22ps_ACOT_NNLO(double const& eta, bool const& xdependent);
    double Regular(double const& x) const;
  };

  class C22g_ACOT_NNLO: public Expression
  {
  public:
    C22g_ACOT_NNLO(double const& eta, bool const& xdependent);
    double Regular(double const& x) const;
    double Local(double const& x)   const;
  };

  //F3
  

  class C32nsNC_ACOT : public Expression
  {
  public:
    C32nsNC_ACOT(int const& nf, double const& eta, bool const& xdependent);
    double Singular(double const& x) const;
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  //FL
  class CL2nsNC_ACOT : public Expression
  {
  public:
    CL2nsNC_ACOT(int const& nf, double const& eta, bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
  private:
    int const _nf;
  };

  class CL2psNC_ACOT : public Expression
  {
  public:
    CL2psNC_ACOT(double const& eta, bool const& xdependent=true);
    double Regular(double const& x)  const;
  };

  class CL2gNC_ACOT : public Expression
  {
  public:
    CL2gNC_ACOT(double const& eta, bool const& xdependent=true);
    double Regular(double const& x)  const;
  };

  // Charged Current

  

  // REMOVE?
  // class C21CCg_general_mass_FRED: public Expression
  // {
  // public:
  //   C21CCg_general_mass_FRED(double const& eta, double const& xi, double const& m1, double const& m2);
  //   double Regular(double const& x)  const;
  // private:
  //   const double _xi;
  //   const double _m1;
  //   const double _m2;
  // };

  

  ////////////////////////////////
  //// ACOT-NNLO stuff
  
}