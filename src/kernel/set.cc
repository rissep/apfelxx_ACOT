//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/set.h"
#include "apfel/distribution.h"
#include "apfel/operator.h"
#include "apfel/tools.h"

namespace apfel {

  //_________________________________________________________________________
  template<class T>
  Set<T>::Set(ConvolutionMap const& map, unordered_map<int,T> const& in):
    _map(map),
    _objects(in)
  {
  }

  //_________________________________________________________________________
  template<class T>
  Set<T>& Set<T>::operator = (Set<T> const& d)
  {
    if(this != &d)
      *this = d;
    return *this;
  }

  //_________________________________________________________________________
  template<class T>
  template<class V> Set<V> Set<T>::operator *= (Set<V> const& d) const
  {
    if (&_map != &d.GetMap())
      throw runtime_exception("Set::operator *=", "Convolution Map does not match (1)");

    unordered_map<int,V> mmap;
    for (auto const& item: _map.GetRules())
      {
        auto o = std::begin(item.second);
        V result = (*o).coefficient * _objects.at((*o).operand) * d.GetObjects().at((*o).object);
        o++;
        for (auto end = std::end(item.second); o != end; o++)
          result += (*o).coefficient * _objects.at((*o).operand) * d.GetObjects().at((*o).object);
        mmap.insert({item.first,result});
      }

    return Set<V>{d.GetMap(),mmap};
  }

  //_________________________________________________________________________
  template<class T>
  Set<T>& Set<T>::operator *= (double const& s)
  {
    for (auto& v: _objects)
      v.second *= s;
    return *this;
  }

  //_________________________________________________________________________
  template<class T>
  Set<T>& Set<T>::operator /= (int const& s)
  {
    const double r = 1.0/static_cast<double>(s);
    for (auto& v: _objects)
      v.second *= r;
    return *this;
  }

  //_________________________________________________________________________
  template<class T>
  Set<T>& Set<T>::operator *= (Set<T> const& d)
  {
    if (&_map != &d.GetMap())
      throw runtime_exception("Set::operator *=", "Convolution Map does not match (2)");

    for (auto& v: _objects)
      v.second *= d.at(v.first);

    return *this;
  }

  //_________________________________________________________________________
  template<class T>
  Set<T>& Set<T>::operator += (Set<T> const& d)
  {
    if (&_map != &d.GetMap())
      throw runtime_exception("Set::operator +=", "Convolution Map does not match");

    for (auto& v: _objects)
      v.second += d.at(v.first);

    return *this;
  }

  template class Set<Distribution>;
  template class Set<Operator>;
  template Set<Distribution> Set<Operator>::operator *= (Set<Distribution> const&) const;
}
