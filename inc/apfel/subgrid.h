//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <vector>

namespace apfel
{
  /**
   * @brief Class for the x-space interpolation SubGrids.
   *
   * Subgrids are the building blocks of the interpolation procedure.
   * This class defines the "SubGrid" object that includes, apart from the
   * grid itself, also the relevant parameters.
   *
   */
  class SubGrid
  {
  public:

    SubGrid() = delete;

    /**
     * @brief Standard internal grid constructor.
     * @param nx number of grid points in x.
     * @param xMin lower edge x of the grid.
     * @param InterDegree interpolation degree.
     */
    SubGrid(int const& nx, double const& xMin, int const& InterDegree);

    /**
     * @brief External grid constructor
     * @param xsg a std::vector with the nodes of the grid
     * @param InterDegree interpolation degree
     */
    SubGrid(std::vector<double> const& xsg, int const& InterDegree);

    /**
     * @brief Check whether SubGrids are equal
     * @param sg the SubGrid to be compared
     * @return true/false
     */
    bool operator == (SubGrid const& sg) const;
    bool operator != (SubGrid const& sg) const;

    // Getters
    int                        nx()          const { return _nx; }          //!< return the number of x points
    int                        InterDegree() const { return _InterDegree; } //!< return the interpolation degree
    bool                       IsExternal()  const { return _IsExternal; }  //!< return external status of sub-grid
    double                     xMin()        const { return _xMin; }        //!< return the minimum node value
    double                     xMax()        const { return _xMax; }        //!< return the maximum node value
    double                     Step()        const { return _Step; }        //!< return the step size of the log grid
    std::vector<double> const& GetGrid()     const { return _xsg; }         //!< return the grid setup.
    std::vector<double> const& GetLogGrid()  const { return _lxsg; }        //!< return the log grid obj.

  private:
    int                 _nx;           //!< Number intervals
    int                 _InterDegree;  //!< Interpolation degree
    bool                _IsExternal;   //!< Is external
    double              _xMin;         //!< Minumim value of x
    double              _xMax;         //!< Maximum value of x (should always be 1)
    double              _Step;         //!< Step pf the logarthmically spaced grid
    std::vector<double> _xsg;          //!< Actual grid
    std::vector<double> _lxsg;         //!< The log of the actual grid

    friend std::ostream& operator<<(std::ostream& os, SubGrid const& dt);
  };

  /**
   * @brief Method which prints SubGrid with cout <<.
   */
  std::ostream& operator<<(std::ostream& os, SubGrid const& sg);
}
