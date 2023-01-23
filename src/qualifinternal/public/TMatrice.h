#ifndef QUALIF_TMATRICE_H
#define QUALIF_TMATRICE_H

#include "Vecteur.h"

namespace Qualif {

template<int DimensionI, int DimensionJ>
class TMatrice
{
private:

  double elements[DimensionI][DimensionJ];
public:

  TMatrice(){};
  TMatrice(Vecteur *v) {
	  for (int i = 0; i <DimensionI; i++)
	      for (int j = 0; j <DimensionJ;j++)
	        elements[i][j] = v[j].GetCoor(i);
  };
  ~TMatrice() {};

  void InitMatrice(int j, Vecteur v)
  {
    for (int i = 0; i <DimensionI; i++)
      elements[i][j] = v.GetCoor(i);
  }

  double Determinant()
  {
    if (DimensionJ == 2)
      {
        Vecteur u(elements[0][0],elements[0][1]);
        Vecteur v(elements[1][0],elements[1][1]);
        return u.Determinant(v);
      }
    else
      {
        Vecteur u(elements[0][0],elements[0][1],elements[0][2]);
        Vecteur v(elements[1][0],elements[1][1],elements[1][2]);
        Vecteur w(elements[2][0],elements[2][1],elements[2][2]);
        return u.Determinant(v,w);
      }
  }

};
}  // end namespace  Qualif
#endif

