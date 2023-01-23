#include "Triangle.h"

namespace Qualif {

//**************************************************************
//             METHODES ET DEFINITIONS DE LA CLASSE 
//                            TRIANGLE
//**************************************************************

//==============================================================
//   CONSTRUCTEUR
//==============================================================
Triangle::Triangle(int DimMail, Vecteur *points) 
{  
  m_DimensionMaillage = DimMail;

  m_Dimension      = 2;
  m_NbSommets      = 3;
  m_NbAretes       = 3;

  Init_Connectivite();
  Init_Sommets(points);
  Init_SommetsIdeaux();
}
Triangle::Triangle(int DimMail)
{
  m_DimensionMaillage = DimMail;

  m_Dimension      = 2;
  m_NbSommets      = 3;
  m_NbAretes       = 3;

  Init_Connectivite();
  Init_SommetsIdeaux();
  m_Sommets   = new Vecteur[m_NbSommets];
}

//==============================================================
//   INITIALISATION DES CONNECTIVITES
//==============================================================
void Triangle::Init_Connectivite()
{
  m_NbVoisins = new int [m_NbSommets];
  m_Voisins   = new int*[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {    
      m_NbVoisins[isommet] = 2;
      m_Voisins  [isommet] = new int[m_NbVoisins[isommet]];
    }
  
  m_Voisins[0][0] = 1; m_Voisins[0][1] = 2;
  m_Voisins[1][0] = 2; m_Voisins[1][1] = 0;
  m_Voisins[2][0] = 0; m_Voisins[2][1] = 1;

  m_SommetsDesAretes = new int*[m_NbAretes];
  for (int iarete = 0; iarete < m_NbAretes; iarete++)
    m_SommetsDesAretes[iarete] = new int[2];

  m_SommetsDesAretes[0][0]=0; m_SommetsDesAretes[0][1]=1;
  m_SommetsDesAretes[1][0]=1; m_SommetsDesAretes[1][1]=2;
  m_SommetsDesAretes[2][0]=2; m_SommetsDesAretes[2][1]=0;
}

//==============================================================
//   INITIALISATION DES SOMMETS
//==============================================================
void Triangle::Init_Sommets(Vecteur *points)
{
  m_Sommets   = new Vecteur[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    m_Sommets[isommet] = points[isommet]; 
}
void Triangle::Init_SommetsIdeaux()
{
  m_SommetsIdeaux = new Vecteur[m_NbSommets];

  m_SommetsIdeaux[0] =   Vecteur(0,0); 
  m_SommetsIdeaux[1] =   Vecteur(1,0); 
  m_SommetsIdeaux[2] =   Vecteur(.5,.5*std::sqrt(3.)); 
}

//==============================================================
//  JACOBIENNE POUR UN TRIANGLE
//==============================================================
Matrice Triangle::Jacobienne(Vecteur *Noeuds)
{
  Matrice J(m_Dimension,m_Dimension);

  if (m_DimensionMaillage==2)
    {
      J.InitMatrice(0,Noeuds[1] - Noeuds[0]);
      J.InitMatrice(1,Noeuds[2] - Noeuds[0]);
    }
  else if (m_DimensionMaillage==3)
    {
      Vecteur *NoeudsNouveauRepere = new Vecteur[3];
      Noeuds[0].NouveauRepere(Noeuds[1],Noeuds[2],
			      NoeudsNouveauRepere);      
      J.InitMatrice(0,NoeudsNouveauRepere[1] - NoeudsNouveauRepere[0]);
      J.InitMatrice(1,NoeudsNouveauRepere[2] - NoeudsNouveauRepere[0]);
      delete [] NoeudsNouveauRepere;
    }
  return J;
}

//==============================================================
//  SURFACE
//==============================================================
double Triangle::VolumeSurface()
{  
  return 0.5 * Jacobienne(m_Sommets).Determinant();
}

//==============================================================
//  MESURE D'ODDY MODIFIEE
//==============================================================
double Triangle::Oddy()
{
  Matrice A(m_Dimension,m_Dimension);
  A = Jacobienne(m_Sommets);
  double Det = A.Determinant();  

  if (Det == 0)
    return Infini;
  else
    {
      double tau;
      Matrice Wm1(m_Dimension,m_Dimension);
      Matrice T  (m_Dimension,m_Dimension);
      
      Wm1    = Jacobienne(m_SommetsIdeaux).inverse(); 
      T      = A*Wm1;
   
      double D1 = std::pow(std::pow(T.Determinant(),4),1./m_Dimension);
      double D2 = std::pow((T.transpose()*T).normeFrobenius(),2);
      double D3 = std::pow(T.normeFrobenius(),4)/m_Dimension;
      return (D2-D3)/D1;
    }
}

//==============================================================
//  CONDITIONNEMENT
//==============================================================
double Triangle::Conditionnement()
{
  Matrice A(m_Dimension,m_Dimension);  
  A = Jacobienne(m_Sommets);
  if (A.Determinant() == 0)
    return Infini;
  else
    return A.normeFrobenius()*(A.inverse()).normeFrobenius()/m_Dimension;
}

//==============================================================
//  SCALEDJACOBIAN
//==============================================================
double Triangle::ScaledJacobian()
{
  Matrice A(m_Dimension,m_Dimension);  
  A = Jacobienne(m_Sommets);
  
  double *L;
  double Lmax;
  L = new double[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {  
      int i0    = m_Voisins[isommet][0];
      int i1    = m_Voisins[isommet][1];
      double l0 =(m_Sommets[isommet] - m_Sommets[i0]).norme2();
      double l1 =(m_Sommets[isommet] - m_Sommets[i1]).norme2();
      L[isommet] = l0*l1;
    }
  Lmax = *PlusGrandElement(L,m_NbSommets);
  delete [] L;
  if (Lmax == 0)
    return Infini;
  else
    return A.Determinant()/ Lmax;
}

//==============================================================
//  KNUPP_SKEW : MESURE D'OBLIQUITE
//==============================================================
double Triangle::Knupp_Skew()
{
  double kappa;
  Matrice A   (m_Dimension,m_Dimension);
  Matrice Qwm1(m_Dimension,m_Dimension);
  Matrice Q   (m_Dimension,m_Dimension);
  Matrice X   (m_Dimension,m_Dimension);

  A  = Jacobienne(m_Sommets); 
  if (A.Determinant() == 0)
    return 0.;
  else
    {
      Qwm1 = Jacobienne(m_SommetsIdeaux).Decompose2D_Q().inverse(); 
      Q    = A.Decompose2D_Q();
      X    = Q*Qwm1;
      kappa =X.normeFrobenius()*(X.inverse().normeFrobenius());
      return (m_Dimension/kappa);
    }
}

//==============================================================
//   KNUPP_SHAPE : MESURE DE FORME
//==============================================================
double Triangle::Knupp_Shape()
{
  Matrice A(m_Dimension,m_Dimension);
  A = Jacobienne(m_Sommets);

  if (A.Determinant() == 0)
    return 0.;
  else
    {  
      double kappa;
      Matrice Wm1(m_Dimension,m_Dimension);
      Matrice T  (m_Dimension,m_Dimension);
      
      Wm1    = Jacobienne(m_SommetsIdeaux).inverse(); 
      T      = A*Wm1;      
      kappa  = T.normeFrobenius()*(T.inverse().normeFrobenius());
      return (m_Dimension/kappa);
    }
}

//==============================================================
//   KNUPP_VOLUME : MESURE DE VOLUME
//==============================================================
double Triangle::Knupp_Volume()
{
  Matrice A(m_Dimension,m_Dimension);
  A = Jacobienne(m_Sommets);

  double Det = A.Determinant();
  if ( Det == 0)
    return 0.;
  else
    {
      Matrice Wm1(m_Dimension,m_Dimension);
      Wm1 = Jacobienne(m_SommetsIdeaux).inverse(); 
      double tau = (A*Wm1).Determinant();
      if (tau > 0)
	return  min(tau,1./tau);
      else
	return  max(tau,1./tau);
      }
}

//==============================================================
//   KNUPP_VOLUMESHAPE : MESURE COMBINEE FORME-VOLUME
//==============================================================
double Triangle::Knupp_VolumeShape()
{
  Matrice A(m_Dimension,m_Dimension);
  A = Jacobienne(m_Sommets);

  double Det = A.Determinant();  
  if (Det == 0)
    return 0.;
  else
    {
      double kappa, tau;
      Matrice Wm1(m_Dimension,m_Dimension);
      Matrice T  (m_Dimension,m_Dimension);
      Wm1    = Jacobienne(m_SommetsIdeaux).inverse(); 
      T      = A*Wm1;
      kappa  = T.normeFrobenius()*(T.inverse().normeFrobenius());
      tau    = T.Determinant();      
      if (tau > 0)
	return  min(tau,1./tau)*m_Dimension/kappa;
      else
	return  max(tau,1./tau)*m_Dimension/kappa;
    }
}

//==============================================================
//   RAPPORT D'ASPECT GAMMA
//==============================================================
double Triangle::AspectRatio_Gamma()
{
  double Surface = VolumeSurface();
  double SommeAretesCarrees = 0.;
  
  for (int iarete = 0; iarete < m_NbAretes; iarete++)
    {
      int i0 = m_SommetsDesAretes[iarete][0];
      int i1 = m_SommetsDesAretes[iarete][1];
      Vecteur arete = m_Sommets[i1] - m_Sommets[i0]; 
      SommeAretesCarrees += arete*arete;
    }
  
  return SommeAretesCarrees/(4*std::sqrt(3.)*Surface);   
}

//==============================================================
//   RAPPORT D'ASPECT Q2
//==============================================================
double Triangle::AspectRatio_Q()
{
  double Surface = VolumeSurface();
  double DemiPerimetre = 0.;
  double RayonInscrit;
  double *L;
  double Lmax;    
  L = new double[m_NbAretes];
  
  for (int iarete = 0; iarete < m_NbAretes; iarete++)
    {
      int i0 = m_SommetsDesAretes[iarete][0];
      int i1 = m_SommetsDesAretes[iarete][1];
      L[iarete]  = (m_Sommets[i1] - m_Sommets[i0]).norme2();
      DemiPerimetre += L[iarete]*0.5;
    }
  Lmax = *PlusGrandElement(L,m_NbAretes);
  delete [] L;
  RayonInscrit = Surface/DemiPerimetre;
  
  return 2*std::sqrt(3.)*RayonInscrit/ Lmax;   
}

//==============================================================
//   ANGLE MINIMUM OU MAXIMUM
//==============================================================
double Triangle::AngleMinMax(const std::string MinMax)
{
  double *angle;
  double angleminmax;  
  Vecteur u,v;
  
  angle = new double[m_NbSommets];

  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {
      int i0 = m_Voisins[isommet][0];
      int i1 = m_Voisins[isommet][1];
      if (m_DimensionMaillage == 2)
	{
	  u  = m_Sommets[i0] - m_Sommets[isommet];
	  v  = m_Sommets[i1] - m_Sommets[isommet];
	}
      else
	{
	  Vecteur *NoeudsNouveauRepere = new Vecteur[3];
	  m_Sommets[isommet].NouveauRepere(m_Sommets[i0],m_Sommets[i1],
			NoeudsNouveauRepere);
	  u  = NoeudsNouveauRepere[1] - NoeudsNouveauRepere[0];
	  v  = NoeudsNouveauRepere[2] - NoeudsNouveauRepere[0];
	  delete [] NoeudsNouveauRepere; 
	}  
      angle[isommet] = u.Angle(v);
     }
  if (MinMax == "Min")
      angleminmax = *PlusPetitElement(angle,m_NbSommets);
  else
      angleminmax = *PlusGrandElement(angle,m_NbSommets);

  delete [] angle;
  return 180./PI*angleminmax; 
}

//==============================================================
//  MINIMUM DES JACOBIENS 
//==============================================================
double Triangle::JacobienMin()
{
  return Jacobienne(m_Sommets).Determinant();
}

}
