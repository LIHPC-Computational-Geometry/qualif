#include "Hexaedre.h"

#include "TMatrice.h"


namespace Qualif {
//**************************************************************
//             METHODES ET DEFINITIONS DE LA CLASSE 
//                            HEXAEDRE
//**************************************************************

//==============================================================
//   CONSTRUCTEUR
//==============================================================
Hexaedre::Hexaedre(Vecteur *points) 
{  
  m_Dimension      =  3;
  m_NbSommets      =  8;
  m_NbAretes       = 12;
  m_NbFaces        =  6;

  Init_Connectivite();
  Init_Sommets(points);
}
Hexaedre::Hexaedre()
{
  m_Dimension      =  3;
  m_NbSommets      =  8;
  m_NbAretes       = 12;
  m_NbFaces        =  6;

  Init_Connectivite();
  m_Sommets   = new Vecteur[m_NbSommets];
}
//==============================================================
//   INITIALISATION DES CONNECTIVITES
//==============================================================
void Hexaedre::Init_Connectivite()
{
  m_NbVoisins = new int [m_NbSommets];
  m_Voisins   = new int*[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {    
      m_NbVoisins[isommet] = 3;
      m_Voisins  [isommet] = new int[m_NbVoisins[isommet]];
    }
  
  m_Voisins[0][0] = 1; m_Voisins[0][1] = 3; m_Voisins[0][2] = 4;
  m_Voisins[1][0] = 2; m_Voisins[1][1] = 0; m_Voisins[1][2] = 5;
  m_Voisins[2][0] = 3; m_Voisins[2][1] = 1; m_Voisins[2][2] = 6;
  m_Voisins[3][0] = 0; m_Voisins[3][1] = 2; m_Voisins[3][2] = 7;
  m_Voisins[4][0] = 7; m_Voisins[4][1] = 5; m_Voisins[4][2] = 0;
  m_Voisins[5][0] = 4; m_Voisins[5][1] = 6; m_Voisins[5][2] = 1;
  m_Voisins[6][0] = 5; m_Voisins[6][1] = 7; m_Voisins[6][2] = 2;
  m_Voisins[7][0] = 6; m_Voisins[7][1] = 4; m_Voisins[7][2] = 3;


  m_SommetsDesAretes = new int*[m_NbAretes];
  for (int iarete = 0; iarete < m_NbAretes; iarete++)
    m_SommetsDesAretes[iarete] = new int[2];

  m_SommetsDesAretes[0][0]=0; m_SommetsDesAretes[0][1]=1;
  m_SommetsDesAretes[1][0]=1; m_SommetsDesAretes[1][1]=2;
  m_SommetsDesAretes[2][0]=2; m_SommetsDesAretes[2][1]=3;
  m_SommetsDesAretes[3][0]=3; m_SommetsDesAretes[3][1]=0;
  m_SommetsDesAretes[4][0]=4; m_SommetsDesAretes[4][1]=5;
  m_SommetsDesAretes[5][0]=5; m_SommetsDesAretes[5][1]=6;
  m_SommetsDesAretes[6][0]=6; m_SommetsDesAretes[6][1]=7;
  m_SommetsDesAretes[7][0]=7; m_SommetsDesAretes[7][1]=4;
  m_SommetsDesAretes[8][0]=0; m_SommetsDesAretes[8][1]=4;
  m_SommetsDesAretes[9][0]=5; m_SommetsDesAretes[9][1]=1;
  m_SommetsDesAretes[10][0]=6; m_SommetsDesAretes[10][1]=2;
  m_SommetsDesAretes[11][0]=7; m_SommetsDesAretes[11][1]=3;

  m_NbSommetsDesFaces= new int [m_NbFaces];
  m_SommetsDesFaces         = new int*[m_NbFaces];
  for (int i = 0; i < m_NbFaces; i++)
    {
      m_NbSommetsDesFaces[i] = 4;
      m_SommetsDesFaces[i]          = new int[m_NbSommetsDesFaces[i]];
    }
      m_SommetsDesFaces[0][0] = 0;  m_SommetsDesFaces[1][0] = 1;  m_SommetsDesFaces[2][0] = 2;  
      m_SommetsDesFaces[0][1] = 1;  m_SommetsDesFaces[1][1] = 2;  m_SommetsDesFaces[2][1] = 6;  
      m_SommetsDesFaces[0][2] = 5;  m_SommetsDesFaces[1][2] = 6;  m_SommetsDesFaces[2][2] = 7;  
      m_SommetsDesFaces[0][3] = 4;  m_SommetsDesFaces[1][3] = 5;  m_SommetsDesFaces[2][3] = 3;  

      m_SommetsDesFaces[3][0] = 0;  m_SommetsDesFaces[4][0] = 0;  m_SommetsDesFaces[5][0] = 4;  
      m_SommetsDesFaces[3][1] = 3;  m_SommetsDesFaces[4][1] = 1;  m_SommetsDesFaces[5][1] = 5;  
      m_SommetsDesFaces[3][2] = 7;  m_SommetsDesFaces[4][2] = 2;  m_SommetsDesFaces[5][2] = 6;  
      m_SommetsDesFaces[3][3] = 4;  m_SommetsDesFaces[4][3] = 3;  m_SommetsDesFaces[5][3] = 7;  
}

//==============================================================
//   INITIALISATION DES SOMMETS
//==============================================================
void Hexaedre::Init_Sommets(Vecteur *points)
{
  m_Sommets   = new Vecteur[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    m_Sommets[isommet] = points[isommet]; 
}

//==============================================================
//  JACOBIENNES AUX SOMMETS 
//==============================================================
Matrice Hexaedre::Jacobienne(int NumSommet,Vecteur *Noeuds)
{
  Matrice J(m_Dimension,m_Dimension);
  switch(NumSommet)
    {
    case 0:
      J.InitMatrice(0,Noeuds[1]-Noeuds[0]);
      J.InitMatrice(1,Noeuds[3]-Noeuds[0]); 
      J.InitMatrice(2,Noeuds[4]-Noeuds[0]); 
      break;	      
    case 1:	      
      J.InitMatrice(0,Noeuds[1]-Noeuds[0]);
      J.InitMatrice(1,Noeuds[2]-Noeuds[1]); 
      J.InitMatrice(2,Noeuds[5]-Noeuds[1]); 
      break;	      
    case 2:	      
      J.InitMatrice(0,Noeuds[2]-Noeuds[3]);
      J.InitMatrice(1,Noeuds[2]-Noeuds[1]); 
      J.InitMatrice(2,Noeuds[6]-Noeuds[2]); 
      break;	      
    case 3:	      
      J.InitMatrice(0,Noeuds[2]-Noeuds[3]);
      J.InitMatrice(1,Noeuds[3]-Noeuds[0]); 
      J.InitMatrice(2,Noeuds[7]-Noeuds[3]); 
      break;
    case 4:	      
      J.InitMatrice(0,Noeuds[5]-Noeuds[4]);
      J.InitMatrice(1,Noeuds[7]-Noeuds[4]); 
      J.InitMatrice(2,Noeuds[4]-Noeuds[0]); 
      break;
    case 5:	      
      J.InitMatrice(0,Noeuds[5]-Noeuds[4]);
      J.InitMatrice(1,Noeuds[6]-Noeuds[5]); 
      J.InitMatrice(2,Noeuds[5]-Noeuds[1]); 
      break;
    case 6:	      
      J.InitMatrice(0,Noeuds[6]-Noeuds[7]);
      J.InitMatrice(1,Noeuds[6]-Noeuds[5]); 
      J.InitMatrice(2,Noeuds[6]-Noeuds[2]); 
      break;
    case 7:	      
      J.InitMatrice(0,Noeuds[6]-Noeuds[7]);
      J.InitMatrice(1,Noeuds[7]-Noeuds[4]); 
      J.InitMatrice(2,Noeuds[7]-Noeuds[3]); 
      break;
    default:
      break;
    }     
  return J;
}
TMatrice<3,3> Hexaedre::TJacobienne(int NumSommet,Vecteur *Noeuds)
{
  TMatrice<3,3> J;
  switch(NumSommet)
    {
    case 0:
      J.InitMatrice(0,Noeuds[1]-Noeuds[0]);
      J.InitMatrice(1,Noeuds[3]-Noeuds[0]);
      J.InitMatrice(2,Noeuds[4]-Noeuds[0]);
      break;
    case 1:
      J.InitMatrice(0,Noeuds[1]-Noeuds[0]);
      J.InitMatrice(1,Noeuds[2]-Noeuds[1]);
      J.InitMatrice(2,Noeuds[5]-Noeuds[1]);
      break;
    case 2:
      J.InitMatrice(0,Noeuds[2]-Noeuds[3]);
      J.InitMatrice(1,Noeuds[2]-Noeuds[1]);
      J.InitMatrice(2,Noeuds[6]-Noeuds[2]);
      break;
    case 3:
      J.InitMatrice(0,Noeuds[2]-Noeuds[3]);
      J.InitMatrice(1,Noeuds[3]-Noeuds[0]);
      J.InitMatrice(2,Noeuds[7]-Noeuds[3]);
      break;
    case 4:
      J.InitMatrice(0,Noeuds[5]-Noeuds[4]);
      J.InitMatrice(1,Noeuds[7]-Noeuds[4]);
      J.InitMatrice(2,Noeuds[4]-Noeuds[0]);
      break;
    case 5:
      J.InitMatrice(0,Noeuds[5]-Noeuds[4]);
      J.InitMatrice(1,Noeuds[6]-Noeuds[5]);
      J.InitMatrice(2,Noeuds[5]-Noeuds[1]);
      break;
    case 6:
      J.InitMatrice(0,Noeuds[6]-Noeuds[7]);
      J.InitMatrice(1,Noeuds[6]-Noeuds[5]);
      J.InitMatrice(2,Noeuds[6]-Noeuds[2]);
      break;
    case 7:
      J.InitMatrice(0,Noeuds[6]-Noeuds[7]);
      J.InitMatrice(1,Noeuds[7]-Noeuds[4]);
      J.InitMatrice(2,Noeuds[7]-Noeuds[3]);
      break;
    default:
      break;
    }
  return J;
}
//==============================================================
//  JACOBIENNE AU CENTRE
//==============================================================
Matrice Hexaedre::JacobienneAuCentre(Vecteur *Noeuds)
{
  Vecteur Nul(0,0,0);
  Matrice J(m_Dimension,m_Dimension);

  J.InitMatrice(0,Nul-Noeuds[0]+Noeuds[1]+Noeuds[2]-Noeuds[3]
		     -Noeuds[4]+Noeuds[5]+Noeuds[6]-Noeuds[7]);
  J.InitMatrice(1,Nul-Noeuds[0]-Noeuds[1]+Noeuds[2]+Noeuds[3]
		     -Noeuds[4]-Noeuds[5]+Noeuds[6]+Noeuds[7]);
  J.InitMatrice(2,Nul-Noeuds[0]-Noeuds[1]-Noeuds[2]-Noeuds[3]
		     +Noeuds[4]+Noeuds[5]+Noeuds[6]+Noeuds[7]);


  return J*0.25;
}
TMatrice<3,3> Hexaedre::TJacobienneAuCentre(Vecteur *Noeuds)
{
  Vecteur Nul(0,0,0);
  TMatrice<3,3> J;

  J.InitMatrice(0,(Nul-Noeuds[0]+Noeuds[1]+Noeuds[2]-Noeuds[3]
		     -Noeuds[4]+Noeuds[5]+Noeuds[6]-Noeuds[7])*0.25);
  J.InitMatrice(1,(Nul-Noeuds[0]-Noeuds[1]+Noeuds[2]+Noeuds[3]
		     -Noeuds[4]-Noeuds[5]+Noeuds[6]+Noeuds[7])*0.25);
  J.InitMatrice(2,(Nul-Noeuds[0]-Noeuds[1]-Noeuds[2]-Noeuds[3]
		     +Noeuds[4]+Noeuds[5]+Noeuds[6]+Noeuds[7])*0.25);

  return J;
}
//==============================================================
//  VOLUME
//==============================================================
double Hexaedre::VolumeSurface()
{  
  double Resu = 0.;
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {
      int i0 = m_Voisins[isommet][0];
      int i1 = m_Voisins[isommet][1];
      int i2 = m_Voisins[isommet][2];
      Vecteur p = m_Sommets[i0] - m_Sommets[isommet];
      Vecteur q = m_Sommets[i1] - m_Sommets[isommet];
      Vecteur r = m_Sommets[i2] - m_Sommets[isommet];
      Resu += p*(q^r)/6;

    }
  Vecteur u0 = m_Sommets[7] - m_Sommets[0];
  Vecteur v0 = m_Sommets[5] - m_Sommets[0];
  Vecteur w0 = m_Sommets[2] - m_Sommets[0];
  Resu += u0*(v0^w0)/6;

  Vecteur u6 = m_Sommets[1] - m_Sommets[6];
  Vecteur v6 = m_Sommets[4] - m_Sommets[6];
  Vecteur w6 = m_Sommets[3] - m_Sommets[6];
  Resu += u6*(v6^w6)/6;
  
		      
  return 0.5*Resu;
}
//==============================================================
//  MESURE D'ODDY MODIFIEE
//==============================================================
double Hexaedre::Oddy()
{
  double *D;
  double Dmax;
  Matrice A(m_Dimension,m_Dimension);
  D = new double[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {
      A = Jacobienne(isommet,m_Sommets);
      double Det = A.Determinant();
      if (Det == 0)
	{
	  delete [] D;
	  return Infini;
	} 
      else
	{
	  double D1 = std::pow(std::pow(Det,4),1./m_Dimension);
	  double D2 = std::pow((A.transpose()*A).normeFrobenius(),2);
	  double D3 = std::pow(A.normeFrobenius(),4)/m_Dimension;
	  D[isommet] =(D2-D3)/D1;
	}
    }
  Dmax = *PlusGrandElement(D,m_NbSommets);
  delete [] D;
  return Dmax;
}

//==============================================================
//  CONDITIONNEMENT
//==============================================================
double Hexaedre::Conditionnement()
{
  double *kappa;
  double kappamax;
  Matrice A(m_Dimension,m_Dimension);  
  kappa = new double[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {  
      A = Jacobienne(isommet,m_Sommets);
      if (A.Determinant() == 0)
	{
	  delete [] kappa;
	  return Infini;
	}
      else
	kappa[isommet] = A.normeFrobenius()*(A.inverse()).normeFrobenius();
    }
  kappamax = *PlusGrandElement(kappa,m_NbSommets);
  delete [] kappa;
  return kappamax/m_Dimension;
}

//==============================================================
//  SCALEDJACOBIAN
//==============================================================
double Hexaedre::ScaledJacobian()
{
  double *scaledJ;
  double scaledJmin;
  scaledJ = new double[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {  
      int i0    = m_Voisins[isommet][0];
      int i1    = m_Voisins[isommet][1];
      int i2    = m_Voisins[isommet][2];
      double l0 =(m_Sommets[isommet] - m_Sommets[i0]).norme2();
      double l1 =(m_Sommets[isommet] - m_Sommets[i1]).norme2();
      double l2 =(m_Sommets[isommet] - m_Sommets[i2]).norme2();
      double l012 = l0*l1*l2;
      double Det = TJacobienne(isommet,m_Sommets).Determinant();
      if (l012 == 0 || Det == 0)
	scaledJ[isommet] = 0.;
      else
	scaledJ[isommet] = Det/l012;
    }
  scaledJmin = *PlusPetitElement(scaledJ,m_NbSommets);
  delete [] scaledJ;

  if(scaledJmin >  1.) scaledJmin =  1.;
  if(scaledJmin < -1.) scaledJmin = -1.;

  return scaledJmin;
}

//==============================================================
//  KNUPP_SKEW : MESURE D'OBLIQUITE
//==============================================================
double Hexaedre::Knupp_Skew()
{
  double *kappa;
  double kappaMax;    
  kappa = new double[m_NbSommets];

  Matrice A(m_Dimension,m_Dimension);
  Matrice Q(m_Dimension,m_Dimension);

  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {
      A = Jacobienne(isommet,m_Sommets);
      if (A.Determinant() == 0)
	{
	  delete [] kappa;
	  return 0.;
	}
      else
	{
	  Q = A.Decompose3D_Q();
	  kappa[isommet] = Q.normeFrobenius()*(Q.inverse().normeFrobenius());
	}
    }
      kappaMax = *PlusGrandElement(kappa,m_NbSommets);
      delete [] kappa;
      return (m_Dimension/ kappaMax);
}

//==============================================================
//   KNUPP_SHAPE : MESURE DE FORME
//==============================================================
double Hexaedre::Knupp_Shape()
{
  double *kappa;
  double kappaMax;      
  Matrice A(m_Dimension,m_Dimension);
  kappa = new double[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {
      A = Jacobienne(isommet,m_Sommets);
      if (A.Determinant() == 0)
	{
	  delete [] kappa;
	  return 0.;
	}
      else
	kappa[isommet] = A.normeFrobenius()*(A.inverse().normeFrobenius());
    }
  kappaMax = *PlusGrandElement(kappa,m_NbSommets);
  delete [] kappa;
  return (m_Dimension/ kappaMax);
}

//==============================================================
//   KNUPP_VOLUME : MESURE DE VOLUME
//==============================================================
double Hexaedre::Knupp_Volume()
{
  double *tau;
  double SommeTau = 0;
  Matrice A(m_Dimension,m_Dimension);
  
  tau = new double[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {
      double Det = TJacobienne(isommet,m_Sommets).Determinant();
      if ( Det == 0)
	tau[isommet] = 0.;
      else
	tau[isommet] = (Det>0)?min(Det,1./Det):max(Det,1./Det);
      SommeTau +=tau[isommet]; 
    }
  delete [] tau;
  return (SommeTau/m_NbSommets) ;
}

//==============================================================
//   KNUPP_VOLUMESHAPE : MESURE COMBINEE FORME-VOLUME
//==============================================================
double Hexaedre::Knupp_VolumeShape()
{
  double *tau;
  double *kappa;
  double kappaMax;
  double SommeTau = 0;
  
  Matrice A(m_Dimension,m_Dimension);
  
  tau   = new double[m_NbSommets];
  kappa = new double[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {
      A = Jacobienne(isommet,m_Sommets);
      double Det = A.Determinant();
      if ( Det == 0)
	{
	  tau  [isommet] = 0.;
	  kappa[isommet] = Infini;
	}
      else
	{
	  tau[isommet] = (Det>0)?min(Det,1./Det):max(Det,1./Det);
	  kappa[isommet] = A.normeFrobenius()*(A.inverse().normeFrobenius()); 
	} 
      SommeTau +=tau[isommet];      
    }
  kappaMax = *PlusGrandElement(kappa,m_NbSommets);
  delete [] tau;
  delete [] kappa;
  return (SommeTau/m_NbSommets)*(m_Dimension/kappaMax) ;
}

//==============================================================
//  MINIMUM DES JACOBIENS AUX SOMMETS ET AU CENTRE 
//==============================================================
double Hexaedre::JacobienMin()
{
  double *j;
  double jmin;
  j = new double[m_NbSommets+1];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    j[isommet] = TJacobienne(isommet,m_Sommets).Determinant();

  j[m_NbSommets] = TJacobienneAuCentre(m_Sommets).Determinant();

  jmin = *PlusPetitElement(j,m_NbSommets+1);
  delete [] j;
  return jmin;
}

//==============================================================
//  JACOBIEN AU CENTRE 
//==============================================================
double Hexaedre::JacobienCentre()
{
  return TJacobienneAuCentre(m_Sommets).Determinant();
}

//==============================================================
//   RAPPORT D'ASPECT
//==============================================================
double Hexaedre::AspectRatioCenter()
{
  Vecteur *CentreFaces;
  double h[3];
  double hmax,hmin;

  CentreFaces = new Vecteur[m_NbFaces] ;
 
  for (int iface = 0; iface < m_NbFaces; iface++)
    {
      CentreFaces[iface] = 0.;  
      for (int isommet = 0; isommet < m_NbSommetsDesFaces[iface]; isommet++)
	  CentreFaces[iface] = CentreFaces[iface] + m_Sommets[m_SommetsDesFaces[iface][isommet]]; 
      CentreFaces[iface] = CentreFaces[iface]/m_NbSommetsDesFaces[iface];
    }
  
  h[0] = (CentreFaces[0] - CentreFaces[2]).norme2();
  h[1] = (CentreFaces[1] - CentreFaces[3]).norme2();
  h[2] = (CentreFaces[4] - CentreFaces[5]).norme2();
  
  hmax = *PlusGrandElement(h,3);
  hmin = *PlusPetitElement(h,3);
  delete [] CentreFaces;
  return hmax/hmin; 
}

//==============================================================
//   MESURE D'OBLIQUITE 
//==============================================================
double Hexaedre::Skew()
{
  Vecteur *CentreFaces, Centre;
  double cosTeta,cosPhi;
  Vecteur a1,a2,a5;
  double  l1,l2,l5;

  CentreFaces = new Vecteur[m_NbFaces] ;
 
  for (int iface = 0; iface < m_NbFaces; iface++)
    {
      CentreFaces[iface] = 0.;  
      for (int isommet = 0; isommet < m_NbSommetsDesFaces[iface]; isommet++)
	  CentreFaces[iface] = CentreFaces[iface] + m_Sommets[m_SommetsDesFaces[iface][isommet]]; 
      CentreFaces[iface] = CentreFaces[iface]/m_NbSommetsDesFaces[iface];
    }  
  Centre       = (m_Sommets[0] + m_Sommets[1] + 
		  m_Sommets[2] + m_Sommets[3] +
		  m_Sommets[4] + m_Sommets[5] +
		  m_Sommets[6] + m_Sommets[7] )/8.;
  
  a1  = CentreFaces[1]-Centre;
  a2  = CentreFaces[2]-Centre;
  a5  = CentreFaces[5]-Centre;

  l1  = a1.norme2();
  l2  = a2.norme2();
  l5  = a5.norme2();

  delete [] CentreFaces;
  if (l1 == 0 || l2 == 0|| l5 == 0)
    return 1;
  else
    {
      cosTeta = std::fabs((a1*a2)/(l1*l2));
      cosPhi  = std::fabs((a2*a5)/(l2*l5));
      return max(cosTeta,cosPhi);
    }
}

//==============================================================
//   TAPERS DE ROBINSON 
//==============================================================
double Hexaedre::Taper()
{
  Vecteur Nul(0,0,0);
  Vecteur C110,C101,C011;

  double tau[3];

  C110 = (Nul - m_Sommets[0] + m_Sommets[1] - m_Sommets[2] + m_Sommets[3]
	      - m_Sommets[4] + m_Sommets[5] - m_Sommets[6] + m_Sommets[7])/8.;

  C101 = (Nul - m_Sommets[0] - m_Sommets[1] + m_Sommets[2] + m_Sommets[3]
	      + m_Sommets[4] + m_Sommets[5] - m_Sommets[6] - m_Sommets[7])/8.;
  
  C011 = (      m_Sommets[0] - m_Sommets[1] - m_Sommets[2] + m_Sommets[3]
	      - m_Sommets[4] + m_Sommets[5] + m_Sommets[6] - m_Sommets[7])/8.;
  
  tau[0] = max(std::fabs(C110.y()),std::fabs(C101.z()));
  tau[1] = max(std::fabs(C110.x()),std::fabs(C011.z()));
  tau[2] = max(std::fabs(C101.x()),std::fabs(C011.y()));

  return *PlusGrandElement(tau,3);
}

//==============================================================
//  ETIREMENT 
//==============================================================
double Hexaedre::Etirement()
{
  double *L;
  double Lmin;
  double Diag[4];

  L = new double [m_NbAretes];
  for (int iarete = 0; iarete < m_NbAretes; iarete++)
    {
      int i0 = m_SommetsDesAretes[iarete][0];
      int i1 = m_SommetsDesAretes[iarete][1];
      L[iarete]  = (m_Sommets[i1] - m_Sommets[i0]).norme2();
    }
  Lmin = *PlusPetitElement(L,m_NbAretes);

  Diag[0] = (m_Sommets[6] - m_Sommets[0]).norme2();
  Diag[1] = (m_Sommets[7] - m_Sommets[1]).norme2();
  Diag[2] = (m_Sommets[4] - m_Sommets[2]).norme2();
  Diag[3] = (m_Sommets[5] - m_Sommets[3]).norme2();

  delete [] L;
  return std::sqrt(3.)*Lmin / *PlusGrandElement(Diag,4);
}


//==============================================================
//  RAPPORT DES DIAGONALES
//==============================================================
double Hexaedre::DiagonalRatio()
{
  double Diag[4];

  Diag[0] = (m_Sommets[6] - m_Sommets[0]).norme2();
  Diag[1] = (m_Sommets[7] - m_Sommets[1]).norme2();
  Diag[2] = (m_Sommets[4] - m_Sommets[2]).norme2();
  Diag[3] = (m_Sommets[5] - m_Sommets[3]).norme2();

  return *PlusPetitElement(Diag,4) / *PlusGrandElement(Diag,4);
}
}

