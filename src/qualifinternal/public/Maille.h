#ifndef DEF_MAILLE_H
#define DEF_MAILLE_H

#include "Vecteur.h"
#include "Matrice.h"

namespace Qualif {

class Maille
{
protected:
  int          m_Dimension;
  int          m_NbSommets, m_NbAretes,m_NbFaces;
  int         *m_NbVoisins;
  int        **m_Voisins;
  int        **m_SommetsDesAretes;
  int         *m_NbSommetsDesFaces;
  int        **m_SommetsDesFaces;
  Vecteur     *m_Sommets;


  Maille();

public:

  virtual ~Maille();

  /// retourne le i-ème noeud de la maille
  virtual Vecteur Sommet (int ind) const;
  virtual void Sommet (int ind, double& x, double& y, double& z) const;

  /// remplace les coordonnées d'un sommet par de nouvelles valeurs
  virtual void Modifier_Sommet(int ind, double x, double y, double z=0.0);

  /// le nombre de sommets
  virtual int NombreSommets() const {return m_NbSommets;}

  /// retourne la dimension de la maille.
  virtual int Dimension () const
  { return m_Dimension; }

  double AppliqueCritere(Critere);

  virtual double VolumeSurface()=0;
  virtual double Paoletti();
  virtual double Oddy()=0;
  virtual double Conditionnement()=0;
  virtual double ScaledJacobian()=0;
  virtual double Knupp_Skew()=0;
  virtual double Knupp_Shape()=0;
  virtual double Knupp_Volume()=0;
  virtual double Knupp_VolumeShape()=0;
  virtual double AspectRatio_Gamma()=0;
  virtual double AspectRatio_Beta()=0;
  virtual double AspectRatio_Q()=0;
  virtual double AngleMinMax(const std::string)=0;
  virtual double JacobienMin()=0;
  virtual double JacobienCentre()=0;
  virtual double AspectRatioCenter()=0;
  virtual double Skew()=0;
  virtual double Taper()=0;
  virtual double Etirement()=0;
  virtual double Warp()=0;
  virtual double WarpBase()=0;
  virtual double DiagonalRatio()=0;
  virtual double Validity() = 0;  // returns either 0.0 or 1.0
};

}
#endif
