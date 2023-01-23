#include "Tetraedre.h"

namespace Qualif {

//**************************************************************
//             METHODES ET DEFINITIONS DE LA CLASSE 
//                            TETRAEDRE
//**************************************************************

//==============================================================
//   CONSTRUCTEUR
//==============================================================
Tetraedre::Tetraedre(Vecteur *points) 
{  
  m_Dimension      = 3;
  m_NbSommets      = 4;
  m_NbAretes       = 6;
  m_NbFaces        = 4;

  Init_Connectivite();
  Init_Sommets(points);
  Init_SommetsIdeaux();
}
Tetraedre::Tetraedre()
{
  m_Dimension      = 3;
  m_NbSommets      = 4;
  m_NbAretes       = 6;
  m_NbFaces        = 4;

  Init_Connectivite();
  Init_SommetsIdeaux();
  m_Sommets   = new Vecteur[m_NbSommets];
}
//==============================================================
//   INITIALISATION DES CONNECTIVITES
//==============================================================
void Tetraedre::Init_Connectivite()
{
  m_NbVoisins = new int [m_NbSommets];
  m_Voisins   = new int*[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {    
      m_NbVoisins[isommet] = 3;
      m_Voisins  [isommet] = new int[m_NbVoisins[isommet]];
    }
  
  m_Voisins[0][0] = 1; m_Voisins[0][1] = 2; m_Voisins[0][2] = 3;
  m_Voisins[1][0] = 2; m_Voisins[1][1] = 0; m_Voisins[1][2] = 3;
  m_Voisins[2][0] = 0; m_Voisins[2][1] = 1; m_Voisins[2][2] = 3;
  m_Voisins[3][0] = 0; m_Voisins[3][1] = 1; m_Voisins[3][2] = 2;


  m_SommetsDesAretes = new int*[m_NbAretes];
  for (int iarete = 0; iarete < m_NbAretes; iarete++)
    m_SommetsDesAretes[iarete] = new int[2];

  m_SommetsDesAretes[0][0] = 0; m_SommetsDesAretes[0][1] = 1;
  m_SommetsDesAretes[1][0] = 1; m_SommetsDesAretes[1][1] = 2; 
  m_SommetsDesAretes[2][0] = 2; m_SommetsDesAretes[2][1] = 0;
  m_SommetsDesAretes[3][0] = 3; m_SommetsDesAretes[3][1] = 0;
  m_SommetsDesAretes[4][0] = 3; m_SommetsDesAretes[4][1] = 1;
  m_SommetsDesAretes[5][0] = 3; m_SommetsDesAretes[5][1] = 2;

  m_NbSommetsDesFaces = new int [m_NbFaces];
  m_SommetsDesFaces             = new int*[m_NbFaces];
  for (int iface = 0; iface < m_NbFaces; iface++)
    {
      m_NbSommetsDesFaces[iface] = 3;
      m_SommetsDesFaces[iface] = new int [m_NbSommetsDesFaces[iface]];
    }
  
  m_SommetsDesFaces[0][0] = 0; m_SommetsDesFaces[0][1] = 1; m_SommetsDesFaces[0][2] = 2;
  m_SommetsDesFaces[1][0] = 0; m_SommetsDesFaces[1][1] = 1; m_SommetsDesFaces[1][2] = 3;
  m_SommetsDesFaces[2][0] = 1; m_SommetsDesFaces[2][1] = 2; m_SommetsDesFaces[2][2] = 3;
  m_SommetsDesFaces[3][0] = 0; m_SommetsDesFaces[3][1] = 3; m_SommetsDesFaces[3][2] = 2;

}

//==============================================================
//   INITIALISATION DES SOMMETS
//==============================================================
void Tetraedre::Init_Sommets(Vecteur *points)
{
  m_Sommets   = new Vecteur[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    m_Sommets[isommet] = points[isommet]; 
}
void Tetraedre::Init_SommetsIdeaux()
{
  m_SommetsIdeaux = new Vecteur[m_NbSommets];

  m_SommetsIdeaux[0] =  Vecteur(0,0,0); 
  m_SommetsIdeaux[1] =  Vecteur(1,0,0); 
  m_SommetsIdeaux[2] =  Vecteur(.5,.5*std::sqrt(3.),0) ; 
  m_SommetsIdeaux[3] =  Vecteur(.5, std::sqrt(3.)/6, std::sqrt(2./3)); 
}

//==============================================================
//  JACOBIENNE POUR UN TETRAEDRE
//==============================================================
Matrice Tetraedre::Jacobienne(Vecteur *Noeuds)
{
  Matrice J(m_Dimension,m_Dimension);
   
  J.InitMatrice(0,Noeuds[1] - Noeuds[0]);
  J.InitMatrice(1,Noeuds[2] - Noeuds[0]);
  J.InitMatrice(2,Noeuds[3] - Noeuds[0]);

  return J;
}

//==============================================================
//  VOLUME
//==============================================================
double Tetraedre::VolumeSurface()
{  
 return Jacobienne(m_Sommets).Determinant()/6;
}

//==============================================================
//  MESURE D'ODDY MODIFIEE
//==============================================================
double Tetraedre::Oddy()
{
  Matrice A(m_Dimension,m_Dimension);
  A = Jacobienne(m_Sommets);
  double Det = A.Determinant();  

  if (Det == 0)
    return Infini;
  else
    {
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
double Tetraedre::Conditionnement()
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
double Tetraedre::ScaledJacobian()
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
      int i2    = m_Voisins[isommet][2];
      double l0 =(m_Sommets[isommet] - m_Sommets[i0]).norme2();
      double l1 =(m_Sommets[isommet] - m_Sommets[i1]).norme2();
      double l2 =(m_Sommets[isommet] - m_Sommets[i2]).norme2();
     L[isommet] = l0*l1*l2;
    }
  Lmax = *PlusGrandElement(L,m_NbSommets);
  delete [] L;
  if (Lmax == 0)
    return 0.;
  else {

	double scaledJmin = (A.Determinant()/ Lmax) * sqrt(2.);
	if(scaledJmin >  1.) scaledJmin =  1.;
	if(scaledJmin < -1.) scaledJmin = -1.;

    return scaledJmin;
  }
}

//==============================================================
//  KNUPP_SKEW : MESURE D'OBLIQUITE
//==============================================================
double Tetraedre::Knupp_Skew()
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
      Qwm1 = Jacobienne(m_SommetsIdeaux).Decompose3D_Q().inverse(); 
      Q    = A.Decompose3D_Q();
      X    = Q*Qwm1;
      kappa =X.normeFrobenius()*(X.inverse().normeFrobenius());
      return (m_Dimension/kappa);
    }
}

//==============================================================
//   KNUPP_SHAPE : MESURE DE FORME
//==============================================================
double Tetraedre::Knupp_Shape()
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
double Tetraedre::Knupp_Volume()
{
  Matrice A(m_Dimension,m_Dimension);
  A = Jacobienne(m_Sommets);

  double Det = A.Determinant();
  if ( Det == 0)
    return 0.;
  else
    {
      Matrice Wm1(m_Dimension,m_Dimension);
      Wm1    = Jacobienne(m_SommetsIdeaux).inverse(); 
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
double Tetraedre::Knupp_VolumeShape()
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
double Tetraedre::AspectRatio_Gamma()
{
  double Volume      = VolumeSurface();
  double VolumeIdeal = Jacobienne(m_SommetsIdeaux).Determinant()/6.;
  double SommeAretesCarrees = 0.;
  
  for (int iarete = 0; iarete < m_NbAretes; iarete++)
    {
      int i0 = m_SommetsDesAretes[iarete][0];
      int i1 = m_SommetsDesAretes[iarete][1];
      Vecteur arete = m_Sommets[i1] - m_Sommets[i0]; 
      SommeAretesCarrees += arete*arete/6.;
    }

  return std::pow(SommeAretesCarrees,3./2.)*VolumeIdeal/ Volume;   
}

//==============================================================
//   RAPPORT D'ASPECT BETA
//==============================================================
double Tetraedre::AspectRatio_Beta()
{
  double Volume           = VolumeSurface();
  double SommeDesSurfaces = 0.;
  double RayonInscrit, RayonCirconscrit;
  
  for (int iface = 0; iface < m_NbFaces; iface++)
    {
      int i0 = m_SommetsDesFaces[iface][0];
      int i1 = m_SommetsDesFaces[iface][1];
      int i2 = m_SommetsDesFaces[iface][2];
      Vecteur *NoeudsNouveauRepere = new Vecteur[3];
      m_Sommets[i0].NouveauRepere(m_Sommets[i1],m_Sommets[i2],
		    NoeudsNouveauRepere);
      Vecteur u = NoeudsNouveauRepere[1] - NoeudsNouveauRepere[0];
      Vecteur v = NoeudsNouveauRepere[2] - NoeudsNouveauRepere[0];
      delete [] NoeudsNouveauRepere;
      SommeDesSurfaces += .5*std::fabs(u.Determinant(v));
    }
  RayonInscrit = 3*Volume/SommeDesSurfaces;

  double a = (m_Sommets[0]-m_Sommets[1]).norme2()*(m_Sommets[2]-m_Sommets[3]).norme2();
  double b = (m_Sommets[0]-m_Sommets[2]).norme2()*(m_Sommets[1]-m_Sommets[3]).norme2();
  double c = (m_Sommets[0]-m_Sommets[3]).norme2()*(m_Sommets[1]-m_Sommets[2]).norme2();
  RayonCirconscrit = std::sqrt((a+b+c)*(a+b-c)*(a+c-b)*(b+c-a))/(24*Volume);

  return RayonCirconscrit/(m_Dimension*RayonInscrit);
 }

//==============================================================
//   RAPPORT D'ASPECT Qmu
//==============================================================
double Tetraedre::AspectRatio_Q()
{
  double Volume           = VolumeSurface();
  double SommeDesSurfaces = 0.;
  double RayonInscrit;
  double *L;
  double Lmax;    
  L = new double[m_NbAretes];
  
  for (int iface = 0; iface < m_NbFaces; iface++)
    {
      int i0 = m_SommetsDesFaces[iface][0];
      int i1 = m_SommetsDesFaces[iface][1];
      int i2 = m_SommetsDesFaces[iface][2];
      Vecteur *NoeudsNouveauRepere = new Vecteur[3];
      m_Sommets[i0].NouveauRepere(m_Sommets[i1],m_Sommets[i2],
		    NoeudsNouveauRepere);
      Vecteur u = NoeudsNouveauRepere[1] - NoeudsNouveauRepere[0];
      Vecteur v = NoeudsNouveauRepere[2] - NoeudsNouveauRepere[0];
      delete [] NoeudsNouveauRepere;
      SommeDesSurfaces += .5*std::fabs(u.Determinant(v));
    }

  RayonInscrit = 3*Volume/SommeDesSurfaces;
  
  for (int iarete = 0; iarete < m_NbAretes; iarete++)
    {
      int i0 = m_SommetsDesAretes[iarete][0];
      int i1 = m_SommetsDesAretes[iarete][1];
      L[iarete]  = (m_Sommets[i1] - m_Sommets[i0]).norme2();
    }
  Lmax = *PlusGrandElement(L,m_NbAretes);
  delete [] L;
  
  return 2*std::sqrt(6.)*RayonInscrit/Lmax;   
}

//==============================================================
//  MINIMUM DES JACOBIENS 
//==============================================================
double Tetraedre::JacobienMin()
{
  return Jacobienne(m_Sommets).Determinant();
}

}

