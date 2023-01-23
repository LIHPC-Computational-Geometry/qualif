#include "Maille.h"
#include <exception>

namespace Qualif {

//**************************************************************
//       METHODES ET DEFINITIONS DE LA CLASSE ABSTRAITE 
//                        MAILLE
//**************************************************************

Maille::Maille()
: m_Dimension(0)
, m_NbSommets(0)
, m_NbAretes(0)
, m_NbFaces(0)
, m_NbVoisins(0)
, m_Voisins(0)
, m_SommetsDesAretes(0)
, m_NbSommetsDesFaces(0)
, m_SommetsDesFaces(0)
, m_Sommets(0)
{
}

Maille::~Maille()
{
	if (m_Sommets)
		delete [] m_Sommets;

	if (m_NbVoisins)
		delete [] m_NbVoisins;
	if (m_Voisins){
		for (int isommet = 0; isommet < m_NbSommets; isommet++)
			delete [] m_Voisins[isommet];
		delete [] m_Voisins;
	}

	if(m_SommetsDesAretes){
		for (int iarete = 0; iarete < m_NbAretes; iarete++)
			delete [] m_SommetsDesAretes[iarete];
		delete [] m_SommetsDesAretes;
	}

	if (m_NbSommetsDesFaces)
		delete [] m_NbSommetsDesFaces;
	if (m_SommetsDesFaces){
		for (int iface = 0; iface < m_NbFaces; iface++)
			delete [] m_SommetsDesFaces[iface];
		delete [] m_SommetsDesFaces;
	}
}

//==============================================================
Vecteur Maille::Sommet(int ind) const
{
	if (ind < 0 || ind >= m_NbSommets)
		throw std::exception();

	return m_Sommets[ind];
}

//==============================================================
void Maille::Sommet(int ind, double& x, double& y, double& z) const
{
	if (ind < 0 || ind >= m_NbSommets)
		throw std::exception();

	x=m_Sommets[ind].x();
	y=m_Sommets[ind].y();
	z=m_Sommets[ind].z();
}

//==============================================================
void Maille::Modifier_Sommet(int ind, double x, double y, double z)
{
	if (ind < 0 || ind >= m_NbSommets)
		throw std::exception();

	m_Sommets[ind].coor[0] = x;
	m_Sommets[ind].coor[1] = y;
	m_Sommets[ind].coor[2] = z;
}
//==============================================================
//    MESURE DE LA QUALITE DE LA MAILLE
//      APPEL AU CRITERE SELECTIONNE 
//==============================================================
double Maille::AppliqueCritere(Critere critere)
{
  switch (critere)
    {
    case VOLUMESURFACE:
      return VolumeSurface();
      break;  
    case PAOLETTI:
      return Paoletti();
      break;  
    case ODDY:
      return Oddy();
      break;  
    case CONDITIONNEMENT:
      return Conditionnement();
      break;  
    case SCALEDJACOBIAN:
      return ScaledJacobian();
      break;  
    case KNUPP_SKEW:
      return Knupp_Skew();
      break;  
    case KNUPP_SHAPE:
      return Knupp_Shape();
      break;  
    case KNUPP_VOLUME:
      return Knupp_Volume();
      break;  
    case KNUPP_VOLUMESHAPE:
      return Knupp_VolumeShape();
      break;  
    case ASPECTRATIO_GAMMA:
      return AspectRatio_Gamma();
      break;  
    case ASPECTRATIO_BETA:
      return AspectRatio_Beta();
      break;  
    case ASPECTRATIO_Q:
      return AspectRatio_Q();
      break;  
    case ANGLEMIN:
      return AngleMinMax("Min");
      break;  
     case ANGLEMAX:
      return AngleMinMax("Max");
      break;  
     case JACOBIENMIN:
      return JacobienMin();
      break;
     case JACOBIENCENTRE:
      return JacobienCentre();
      break;
    case ASPECTRATIOCENTER:
      return AspectRatioCenter();
      break;  
    case SKEW:
      return Skew();
      break;  
    case TAPER:
      return Taper();
      break;  
    case ETIREMENT:
      return Etirement();
      break;  
    case WARP:
      return Warp();
      break;  
    case WARPBASE:
      return WarpBase();
      break; 
    case DIAGONALRATIO:
      return DiagonalRatio();
      break;
    case VALIDITY:
      return Validity();
      break;
    default:
      std::cout << "*** CRITERE NON DEFINI POUR CE TYPE DE MAILLE  ***" << std::endl;
      return Not_A_Number;
      break;
    }
}

//==============================================================
//  MESURE DE PAOLETTI
//==============================================================
double Maille::Paoletti()
{  
  double D; 
  Vecteur *E;
  E = new Vecteur[m_NbAretes];
  for (int iarete = 0; iarete < m_NbAretes; iarete++)
    { 
      int s0 = m_SommetsDesAretes[iarete][0];
      int s1 = m_SommetsDesAretes[iarete][1];
      E[iarete] = m_Sommets[s1] - m_Sommets[s0];
    }

  Matrice J(m_Dimension,m_NbAretes,E);
  Matrice A(m_Dimension,m_Dimension);
  A = J*(J.transpose());
  delete [] E;
  if (A.Determinant() == 0)
    return 0.;
  else 
    return m_Dimension/(A.normeFrobenius()*A.inverse().normeFrobenius());      
}

}

