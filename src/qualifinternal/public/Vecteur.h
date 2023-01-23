#ifndef DEF_VECTEUR_H
#define DEF_VECTEUR_H

#include "IQualif.h"


namespace Qualif {

class Maille;

class Vecteur
{
	friend class Maille;

private:
  double coor[DIM];
public:
  Vecteur(double x = 0, double y = 0, double z = 0);
  Vecteur(const Vecteur& v);
  Vecteur& operator = (const Vecteur& v);
  double x();
  double y();
  double z();
  double GetCoor(int i);
  void   NouveauRepere(Vecteur,Vecteur,Vecteur *);
  inline double  operator*(Vecteur);
  inline Vecteur operator^(Vecteur);
  inline Vecteur operator*(double);
  inline Vecteur operator/(double);
  inline Vecteur operator+(Vecteur);
  inline Vecteur operator-(Vecteur);
  inline double  norme2();
  inline double  Determinant(Vecteur);
  inline double  Determinant(Vecteur,Vecteur);
  inline double  Angle(Vecteur);
};

//==============================================================
//   OPERATEUR * : PRODUIT SCALAIRE
//==============================================================
inline double Vecteur::operator*(Vecteur v)
{
  return (coor[0]*v.coor[0]+coor[1]*v.coor[1]+coor[2]*v.coor[2]);
}

//==============================================================
//   OPERATEUR ^ : PRODUIT VECTORIEL
//==============================================================
inline Vecteur  Vecteur::operator^(Vecteur v)
{
  return Vecteur(coor[1]*v.coor[2] - coor[2]*v.coor[1],  
		 coor[2]*v.coor[0] - coor[0]*v.coor[2],
		 coor[0]*v.coor[1] - v.coor[0]*coor[1]);
}

//==============================================================
//   OPERATEUR * : PRODUIT REEL VECTEUR
//==============================================================
inline Vecteur  Vecteur::operator*(double a)
{
  return Vecteur(a*coor[0],a*coor[1],a*coor[2]);
}

//==============================================================
//   OPERATEUR / : DIVISION VECTEUR REEL
//==============================================================
inline Vecteur  Vecteur::operator/(double a)
{
  return Vecteur(coor[0]/a,coor[1]/a,coor[2]/a);
}

//==============================================================
//   OPERATEUR + 
//==============================================================
inline Vecteur  Vecteur::operator+(Vecteur v)
{
  return Vecteur(coor[0] + v.coor[0],
		 coor[1] + v.coor[1],
		 coor[2] + v.coor[2]);
}
//==============================================================
//   OPERATEUR -
//==============================================================
inline Vecteur  Vecteur::operator-(Vecteur v)
{
  return Vecteur(coor[0] - v.coor[0],
		 coor[1] - v.coor[1],
		 coor[2] - v.coor[2]);
}
//==============================================================
//   NORME2
//==============================================================
inline double Vecteur::norme2()
{
  return std::sqrt((*this) * (*this));
}
 
//==============================================================
//   DETERMINANT 
//==============================================================
inline double Vecteur::Determinant(Vecteur v)
{
  return(coor[0] * v.coor[1] - coor[1] * v.coor[0]);
}
inline double Vecteur::Determinant(Vecteur v, Vecteur w)
{
  return((*this)*(v^w));
}

//==============================================================
//  ANGLE ENTRE DEUX VECTEURS
//  Angle compris entre 0 et 2*Pi 
//==============================================================
inline double Vecteur::Angle(Vecteur v)
{ 
  double ProdScalaire  =  (*this)*v;
  double Surface       = ((*this)^v).z();
  double ProdNormes    = (this->norme2())*v.norme2();

  bool CosPositif = (ProdScalaire > 0.);
  bool SinPositif = (Surface      > 0.);

  if (Surface      == 0.)                return  0.;
  if (ProdScalaire == 0. &&  SinPositif) return  PI*0.5;
  if (ProdScalaire == 0. && !SinPositif) return -PI*0.5;

  if ( CosPositif &&  SinPositif) return      std::asin(std::fabs(Surface)/ProdNormes);
  if (!CosPositif &&  SinPositif) return   PI-std::asin(std::fabs(Surface)/ProdNormes);
  if (!CosPositif && !SinPositif) return  -PI+std::asin(std::fabs(Surface)/ProdNormes);
  
  // assert ( CosPositif && !SinPositif);
  return     -std::asin(std::fabs(Surface)/ProdNormes);
}

}
#endif

 
