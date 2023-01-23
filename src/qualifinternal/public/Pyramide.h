#ifndef DEF_PYRAMIDE_H
#define DEF_PYRAMIDE_H

#include "Maille.h"

namespace Qualif {

class Pyramide : public Maille
{
private :
  Vecteur *m_SommetsIdeaux;
public :
  Pyramide(Vecteur *);
  Pyramide();
  void Init_Connectivite();
  void Init_Sommets(Vecteur *);
  void Init_SommetsIdeaux();
  Matrice Jacobienne(int, Vecteur *);  

  virtual double VolumeSurface();
  virtual double Oddy();
  virtual double Conditionnement();
  virtual double ScaledJacobian(){return Not_A_Number;}
  virtual double Knupp_Skew();
  virtual double Knupp_Shape();
  virtual double Knupp_Volume();
  virtual double Knupp_VolumeShape();
  virtual double AspectRatio_Gamma(){return Not_A_Number;}
  virtual double AspectRatio_Beta(){return Not_A_Number;}
  virtual double AspectRatio_Q(){return Not_A_Number;}
  virtual double AngleMinMax(const std::string){return Not_A_Number;}
  virtual double JacobienMin();
  virtual double JacobienCentre(){return Not_A_Number;}
  virtual double AspectRatioCenter(){return Not_A_Number;}
  virtual double Skew(){return Not_A_Number;}
  virtual double Taper(){return Not_A_Number;}
  virtual double Etirement(){return Not_A_Number;}
  virtual double Warp(){return Not_A_Number;}
  virtual double WarpBase();
  virtual double DiagonalRatio(){return Not_A_Number;}
  virtual double Validity() { return Not_A_Number; }
};

}

#endif
