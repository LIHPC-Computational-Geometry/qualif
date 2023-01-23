#ifndef DEF_HEXAEDRE_H
#define DEF_HEXAEDRE_H

#include "Maille.h"

#include "TMatrice.h"

namespace Qualif {

class Hexaedre : public Maille
{
public :
  Hexaedre(Vecteur *);
  Hexaedre();
  void Init_Connectivite();
  void Init_Sommets(Vecteur *);
  Matrice Jacobienne(int, Vecteur *);  
  TMatrice<3,3> TJacobienne(int, Vecteur *);
  Matrice JacobienneAuCentre(Vecteur *);
  TMatrice<3,3> TJacobienneAuCentre(Vecteur *);


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
  virtual double AngleMinMax(const std::string){return Not_A_Number;};
  virtual double JacobienMin();
  virtual double JacobienCentre();
  virtual double AspectRatioCenter();
  virtual double Skew();
  virtual double Taper();
  virtual double Etirement();
  virtual double Warp(){return Not_A_Number;};
  virtual double WarpBase(){return Not_A_Number;};
  virtual double DiagonalRatio();
  virtual double Validity();
};

}
#endif
