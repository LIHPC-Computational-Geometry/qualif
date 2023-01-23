#ifndef DEF_TRIANGLE_H
#define DEF_TRIANGLE_H

#include "Maille.h"

namespace Qualif {

class Triangle : public Maille
{
private :
  int m_DimensionMaillage;
  Vecteur *m_SommetsIdeaux;
public :
  Triangle(int,Vecteur *);
  Triangle(int);
  void Init_Connectivite();
  void Init_Sommets(Vecteur *);
  void Init_SommetsIdeaux();
  virtual Matrice Jacobienne(Vecteur *);  

  virtual double VolumeSurface();
  virtual double Oddy();
  virtual double Conditionnement();
  virtual double ScaledJacobian();
  virtual double Knupp_Skew();
  virtual double Knupp_Shape();
  virtual double Knupp_Volume();
  virtual double Knupp_VolumeShape();
  virtual double AspectRatio_Gamma();
  virtual double AspectRatio_Beta(){return Not_A_Number;}
  virtual double AspectRatio_Q();
  virtual double AngleMinMax(const std::string);
  virtual double JacobienMin();
  virtual double JacobienCentre(){return Not_A_Number;}
  virtual double AspectRatioCenter(){return Not_A_Number;}
  virtual double Skew(){return Not_A_Number;}
  virtual double Taper(){return Not_A_Number;}
  virtual double Etirement(){return Not_A_Number;}
  virtual double Warp(){return Not_A_Number;}
  virtual double WarpBase(){return Not_A_Number;}
  virtual double DiagonalRatio(){return Not_A_Number;}
  virtual double Validity() {return (JacobienMin() > 0.0) ? 1.0 : 0.0; }
};

}
#endif
