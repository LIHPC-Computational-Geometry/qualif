#ifndef DEF_MATRICE_H
#define DEF_MATRICE_H

#include "Vecteur.h"

namespace Qualif {

class Matrice
{
private:
  int DimensionI;
  int DimensionJ;
  double **elements;
public:
  Matrice(int dimI, int dimJ);
  Matrice(int dimI, int dimJ, Vecteur *v);
  ~Matrice();
  int     DimI();
  int     DimJ();
  Vecteur Colonne(int);
  double  Element(int , int);
  void    InitMatrice(int , Vecteur );
  void    InitMatrice(int , int , double );

  inline Matrice(const Matrice &);
  inline Matrice operator=(Matrice );
  inline Matrice operator*(Matrice );
  inline Matrice operator*(double );
  inline Matrice operator+(Matrice );
  inline Matrice operator-(Matrice );
  inline double  Determinant();
  inline Matrice transpose();
  inline Matrice inverse();
  inline double  normeFrobenius();
  inline Matrice Decompose2D_Q();
  inline Matrice Decompose3D_Q();
};

//==============================================================
//   FONCTIONS INLINE
//==============================================================

//==============================================================
//   OPERATEUR DE COPIE 
//==============================================================
inline  Matrice::Matrice(const Matrice &B)
{
  DimensionI= B.DimensionI; 
  DimensionJ= B.DimensionJ;

  elements = new double* [DimensionI];
  for (int i = 0; i<DimensionI; i++)
    elements[i] = new double[DimensionJ];

  for (int i = 0; i < DimensionI ; i++)
    for (int j = 0; j <  DimensionJ; j++)
      elements[i][j] = B.elements[i][j];
 }
//==============================================================
//   OPERATEUR = 
//==============================================================
inline Matrice Matrice::operator=(Matrice B)
{
  DimensionI= B.DimensionI; 
  DimensionJ= B.DimensionJ;


  for (int i = 0; i < DimensionI ; i++)
    for (int j = 0; j <  DimensionJ; j++)
      elements[i][j] = B.elements[i][j];
  return *this;
}
//==============================================================
//   OPERATEUR * : PRODUIT DE MATRICE
//==============================================================
inline Matrice Matrice::operator*(Matrice B)
{
  int dimI = DimensionI; 
  int dimJ = B.DimensionJ;
  double somme;
  Matrice C(dimI,dimJ);

  for (int i = 0; i <dimI; i++)
    for (int j = 0; j <dimJ; j++)
      {
	somme = 0.;
	for (int k = 0; k < DimensionJ; k++)
	  somme += elements[i][k]*B.elements[k][j];
	C.InitMatrice(i,j,somme);
      }      
  return C;
}
//==============================================================
//   OPERATEUR * : PRODUIT MATRICE-REEL
//==============================================================
inline Matrice Matrice::operator*(double a)
{
  int dimI = DimensionI; 
  int dimJ = DimensionJ;
  double somme;
  Matrice C(dimI,dimJ);
 
  for (int i = 0; i <dimI; i++)
    for (int j = 0; j <dimJ; j++)
      C.InitMatrice(i,j,a*elements[i][j]); 
  
  return C;
}//==============================================================
//   OPERATEUR + 
//==============================================================
inline Matrice Matrice::operator+(Matrice B)
{
  int dimI = DimensionI; 
  int dimJ = B.DimensionJ; 
  Matrice C(dimI,dimJ);

  for (int i = 0; i <dimI; i++)
    for (int j = 0; j <dimJ; j++)
      C.InitMatrice(i,j,elements[i][j]+B.elements[i][j]); 
  return C;
}
//==============================================================
//   OPERATEUR -
//==============================================================
inline Matrice Matrice::operator-(Matrice B)
{
  int dimI = DimensionI; 
  int dimJ = B.DimensionJ; 
  Matrice C(dimI,dimJ);

  for (int i = 0; i <dimI; i++)
    for (int j = 0; j <dimJ; j++)
      C.InitMatrice(i,j,elements[i][j]-B.elements[i][j]); 
  return C;
}
//==============================================================
//   CALCUL DU DETERMINANT (MAXIMUM MATRICE 3*3!!!)
//==============================================================
inline double  Matrice::Determinant()
{
  int dimI = DimensionI; 
  int dimJ = DimensionJ;
  if (dimJ == 2)
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
//==============================================================
//   CALCUL DE LA TRANSPOSEE
//==============================================================
inline Matrice Matrice::transpose()
{
  int dimI = DimensionI; 
  int dimJ = DimensionJ; 
  Matrice tM(dimJ,dimI);
  for (int i = 0; i < dimI; i++)
      for (int j = 0; j < dimJ; j++)
	tM.InitMatrice(j,i,elements[i][j]);
  return tM;
}
//==============================================================
//   CALCUL DE L'INVERSE (MAXIMUM MATRICE 3*3!!!)
//==============================================================
inline Matrice Matrice::inverse()
{
  int dimI = DimensionI; 
  int dimJ = DimensionJ; 
  if(dimI != dimJ) 
    {
      std::cout << "MATRICE NON CARREE : INVERSE IMPOSSIBLE" << std::endl;
      return Matrice (dimI,dimJ);
    }
  
  double delta = this->Determinant();
  if (dimI == 2)
      {
	Vecteur v[2] = {Vecteur( elements[1][1]/delta,-elements[1][0]/delta),
			Vecteur(-elements[0][1]/delta, elements[0][0]/delta)};
	return Matrice(dimI,dimJ,v);
      }
  else if (dimI == 3)
    {
      Vecteur v[3] = {Vecteur(elements[1][1]*elements[2][2] - elements[2][1]*elements[1][2],
			      elements[2][1]*elements[0][2] - elements[0][1]*elements[2][2],
			      elements[0][1]*elements[1][2] - elements[1][1]*elements[0][2]),
		      Vecteur(elements[2][0]*elements[1][2] - elements[1][0]*elements[2][2],
			      elements[0][0]*elements[2][2] - elements[2][0]*elements[0][2],
			      elements[1][0]*elements[0][2] - elements[0][0]*elements[1][2]),
                      Vecteur(elements[1][0]*elements[2][1] - elements[2][0]*elements[1][1],
			      elements[2][0]*elements[0][1] - elements[0][0]*elements[2][1],
			      elements[0][0]*elements[1][1] - elements[1][0]*elements[0][1]),
			};
      for (int i = 0; i < dimI; i++) 
	v[i] = v[i]/delta;
      Matrice Cofacteurs(dimI,dimJ,v);
      return Cofacteurs.transpose();
    }
  return Matrice (dimI,dimJ);
}
//==============================================================
//   NORME MATRICIELLE DE FROBENIUS
//==============================================================
inline double Matrice::normeFrobenius()
{
  int dimI = DimensionI; 
  int dimJ = DimensionJ; 
  double x=0.;
  for (int i = 0; i < dimI; i++)
    for (int j = 0; j < dimJ; j++)
      x += elements[i][j] * elements[i][j];  
  return std::sqrt(x);
}

//==============================================================
// POUR LES CRITERES KNUPP
//==============================================================
//==============================================================
// CALCUL DE LA MATRICE Q (resp. Qw) 
// POUR LES ELEMENTS 2D
//=============================================================
inline Matrice  Matrice::Decompose2D_Q()
{
  int Dim = DimensionI;
  Matrice Lambda(Dim,Dim);
  Lambda = this->transpose() * (*this);
  double denom = std::sqrt( Lambda.Element(0,0)*Lambda.Element(1,1));

  Vecteur u[2] = {Vecteur(1,0), 
		  Vecteur(Lambda.Element(0,1)/denom,this->Determinant()/denom)};

  return Matrice(Dim,Dim,u);
}

//==============================================================
// CALCUL DE LA MATRICE Q (resp. Qw) 
// POUR LES ELEMENTS 3D
//==============================================================
inline Matrice  Matrice::Decompose3D_Q()
{
  int Dim = DimensionI;
  Matrice Lambda(Dim,Dim);
  Lambda = this->transpose()* (*this);

  double x      = (Colonne(0)^Colonne(1)).norme2();
  double denom01 = std::sqrt( Lambda.Element(0,0)*Lambda.Element(1,1));
  double denom02 = std::sqrt( Lambda.Element(0,0)*Lambda.Element(2,2));
  double denom22 = std::sqrt( Lambda.Element(2,2));
  double num    = Lambda.Element(0,0)*Lambda.Element(1,2) -
                  Lambda.Element(0,1)*Lambda.Element(0,2);

  Vecteur u[3] = {Vecteur(1,0,0), 
		  Vecteur(Lambda.Element(0,1)/denom01, x/denom01,0),
		  Vecteur(Lambda.Element(0,2)/denom02,
			  num/(denom02*x),
			  this->Determinant()/(denom22*x))};

  return Matrice(Dim,Dim,u);
}
}
#endif
   
