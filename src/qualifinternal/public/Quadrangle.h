#ifndef DEF_QUADRANGLE_H
#define DEF_QUADRANGLE_H

#include "Maille.h"

namespace Qualif {

class Quadrangle : public Maille
{
private :
  int m_DimensionMaillage;
public :
  Quadrangle(int,Vecteur *);
  Quadrangle(int);
  void Init_Connectivite();
  void Init_Sommets(Vecteur *);
  Matrice Jacobienne(int, Vecteur *);  
  Matrice JacobienneAuCentre(Vecteur *);

  virtual double VolumeSurface();
  virtual double Oddy();
  virtual double Conditionnement();
  virtual double ScaledJacobian();
  virtual double Knupp_Skew();
  virtual double Knupp_Shape();
  virtual double Knupp_Volume();
  virtual double Knupp_VolumeShape();
  virtual double AspectRatio_Gamma(){return Not_A_Number;};
  virtual double AspectRatio_Beta(){return Not_A_Number;};
  virtual double AspectRatio_Q(){return Not_A_Number;};
  virtual double AngleMinMax(const std::string);
  virtual double JacobienMin();
  virtual double JacobienCentre();
  virtual double AspectRatioCenter();
  virtual double Skew();
  virtual double Taper();
  virtual double Etirement();
  virtual double Warp();
  virtual double WarpBase(){return Not_A_Number;};
  virtual double DiagonalRatio(){return Not_A_Number;};
  virtual double Validity() {return (JacobienMin() >= 0.0) ? 1.0 : 0.0; }
};

}
#endif
