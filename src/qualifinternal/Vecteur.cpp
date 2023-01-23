#include "Vecteur.h"

namespace Qualif {

//**************************************************************
//             METHODES ET DEFINITIONS DE LA CLASSE 
//                            VECTEUR
//**************************************************************

//==============================================================
//   Constructeur
//==============================================================
Vecteur::Vecteur(double x, double y, double z)
{
  coor[0] = x;
  coor[1] = y;
  coor[2] = z;
}

//==============================================================
//   Constructeur de copie
//==============================================================
Vecteur::Vecteur(const Vecteur& v)
{
  for (int i = 0; i < 3; i++)
     coor[i] = v.coor[i];
}

//==============================================================
//   Opérateur d'affectation
//==============================================================
Vecteur& Vecteur::operator =(const Vecteur& v)
{
  if (&v != this)
  {
    for (int i = 0; i < 3; i++)
       coor[i] = v.coor[i];
  }

  return *this;
}

//==============================================================
//   RENVOI LES COORDONNEES
//==============================================================
double Vecteur::x(){return coor[0];}
double Vecteur::y(){return coor[1];}
double Vecteur::z(){return coor[2];}
double Vecteur::GetCoor(int i){return coor[i];}

//==============================================================
// CHANGEMENT DE COORDONNEES (TRIANGLE ET QUADRANGLE 3D...)
//
// Si le maillage est de dimension 3
// on passe toutes les coordonnees dans le plan du triangle
// Pour orienter le nouveau triangle on regarde si le triedre
// defini par les 3 sommets et un 4eme point situe en moins 
// l'infini(ici : (-1e4,-1e4,-1e4)) est direct.
//==============================================================
void Vecteur::NouveauRepere(
			  Vecteur Sommet1,
			  Vecteur Sommet2,
			  Vecteur *SommetsNouveauRepere)
{
  double L0 = (Sommet1- *this   ).norme2();
  double L1 = (Sommet2-Sommet1).norme2();
  double L2 = (*this-Sommet2   ).norme2();
  double alpha = (L0*L0 - L1*L1 + L2*L2)/ (2*L0);
  double beta  = std::sqrt(L2*L2 - alpha*alpha);
  
  Vecteur MoinsInfini(-1e4,-1e4,-1e4);

  SommetsNouveauRepere[0] = Vecteur(0,0);
  SommetsNouveauRepere[1] = Vecteur(L0,0);
  if ( ((*this-MoinsInfini)^(Sommet1-MoinsInfini))
       *(Sommet2-MoinsInfini) >= 0)
    SommetsNouveauRepere[2] = Vecteur(alpha,beta);
  else
    SommetsNouveauRepere[2] = Vecteur(alpha,-beta);
}

}
