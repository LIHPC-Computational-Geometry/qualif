#ifndef DEF_CONSTQUALIF_H 
#define DEF_CONSTQUALIF_H

#include <cmath>	// asin

namespace Qualif {

//==============================================================
//    NOMBRE DE COORDONNEES DES VECTEURS
//==============================================================
const int DIM = 3;

//==============================================================
//    PI
//==============================================================
const double PI = 2*std::asin(1.);

extern const double Not_A_Number, Infini;

//==============================================================
//     LES CRITERES DISPONIBLES DANS LA BIBLIOTHEQUE
//==============================================================
enum Critere
{
  VOLUMESURFACE,
  PAOLETTI,
  ODDY,
  CONDITIONNEMENT,
  SCALEDJACOBIAN,
  KNUPP_SKEW,
  KNUPP_SHAPE,
  KNUPP_VOLUME,
  KNUPP_VOLUMESHAPE,
  ASPECTRATIO_GAMMA,
  ASPECTRATIO_BETA,
  ASPECTRATIO_Q,
  ANGLEMIN,
  ANGLEMAX,
  JACOBIENMIN,
  JACOBIENCENTRE,
  ASPECTRATIOCENTER,
  SKEW,
  TAPER,
  ETIREMENT,
  WARP,
  WARPBASE,
  DIAGONALRATIO,
  VALIDITY,
  FIN
};
const std::string CRITERESTR[] =
{
  "VolumeSurface",
  "Paoletti",
  "Oddy",
  "Conditionnement",
  "ScaledJacobian",
  "Knupp_Skew",
  "Knupp_Shape",
  "Knupp_Volume",
  "Knupp_VolumeShape",
  "AspectRatio_Gamma",
  "AspectRatio_Beta",
  "AspectRatio_Q",
  "AngleMin",
  "AngleMax",
  "JacobienMin",
  "JacobienCentre",
  "AspectRatioCenter",
  "Skew",
  "Taper",
  "Etirement",
  "Warp",
  "WarpBase",
  "DiagonalRatio",
  "Validity"
};

}

#endif
