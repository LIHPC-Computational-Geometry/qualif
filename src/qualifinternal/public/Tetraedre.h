#ifndef DEF_TETRAEDRE_H
#define DEF_TETRAEDRE_H

#include "Maille.h"

namespace Qualif {
class Tetraedre : public Maille
{
private :
  Vecteur *m_SommetsIdeaux;
public :
  Tetraedre(Vecteur *);
  Tetraedre();
  void Init_Connectivite();
  void Init_Sommets(Vecteur *);
  void Init_SommetsIdeaux();
  Matrice Jacobienne(Vecteur *);  


  virtual double VolumeSurface();
  virtual double Oddy();
  virtual double Conditionnement();
  virtual double ScaledJacobian();
  virtual double Knupp_Skew();
  virtual double Knupp_Shape();
  virtual double Knupp_Volume();
  virtual double Knupp_VolumeShape();
  virtual double AspectRatio_Gamma();
  virtual double AspectRatio_Beta();
  virtual double AspectRatio_Q();
  virtual double AngleMinMax(const std::string){return Not_A_Number;}
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
