#include "Matrice.h"


//**************************************************************
//             METHODES ET DEFINITIONS DE LA CLASSE 
//                            MATRICE
//**************************************************************

namespace Qualif {

//==============================================================
//   CONSTRUCTEURS
//==============================================================
Matrice::Matrice(int dimI, int dimJ)
{
  DimensionI = dimI;
  DimensionJ = dimJ;
  elements = new double* [DimensionI];
  for (int i = 0; i<DimensionI; i++)
    elements[i] = new double[DimensionJ];

  for (int i = 0; i <DimensionI; i++)
    for (int j = 0; j <DimensionJ; j++)
      elements[i][j] = 0.;
}
Matrice::Matrice(int dimI, int dimJ, Vecteur *v)
{
  DimensionI = dimI;
  DimensionJ = dimJ;
  elements = new double*[DimensionI];
  for (int i = 0; i<DimensionI; i++)
    elements[i] = new double[DimensionJ];

  for (int i = 0; i <DimensionI; i++)
    for (int j = 0; j <DimensionJ;j++)
      elements[i][j] = v[j].GetCoor(i);
}

//==============================================================
//   DESTRUCTEUR
//==============================================================
Matrice::~Matrice() 
{
  for (int i = 0; i<DimensionI; i++)
    delete[] elements[i] ;

  delete [] elements;
}
//==============================================================
//   RENVOI DES DIMENSIONS
//==============================================================
int Matrice::DimI(){return DimensionI;}
int Matrice::DimJ(){return DimensionJ;}

//==============================================================
//   RENVOI UN ELEMENT
//==============================================================
double Matrice::Element(int i, int j)
{
  return elements[i][j];
}

//==============================================================
//   RENVOI UN VECTEUR CORRESPONDANT A LA JEME COLONNE
//==============================================================
Vecteur Matrice::Colonne(int j)
{
  if (DimensionI == 2)
    return Vecteur(elements[0][j],elements[1][j]);
  else if (DimensionI == 3)
    return Vecteur(elements[0][j],elements[1][j],elements[2][j]);    
}

//==============================================================
//   INITIALISE UN ELEMENT
//==============================================================
void Matrice::InitMatrice(int i, int j, double val)
{
  elements[i][j] = val;
}
//==============================================================
//   INITIALISE UNE COLONNE
//==============================================================
void Matrice::InitMatrice(int j, Vecteur v)
{
  for (int i = 0; i <DimensionI; i++)
    elements[i][j] = v.GetCoor(i);
}


}

