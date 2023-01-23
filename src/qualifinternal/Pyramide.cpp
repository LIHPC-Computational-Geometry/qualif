#include "Pyramide.h"
#include "Quadrangle.h"

namespace Qualif {

//**************************************************************
//             METHODES ET DEFINITIONS DE LA CLASSE 
//                            PYRAMIDE
//**************************************************************

//==============================================================
//   CONSTRUCTEUR
//==============================================================
Pyramide::Pyramide(Vecteur *points) 
{  
  m_Dimension      = 3;
  m_NbSommets      = 5;
  m_NbAretes       = 8;
  m_NbFaces        = 5;

  Init_Connectivite();
  Init_Sommets(points);
  Init_SommetsIdeaux();
}
Pyramide::Pyramide()
{
  m_Dimension      = 3;
  m_NbSommets      = 5;
  m_NbAretes       = 8;
  m_NbFaces        = 5;

  Init_Connectivite();
  Init_SommetsIdeaux();
  m_Sommets   = new Vecteur[m_NbSommets];
}
//==============================================================
//   INITIALISATION DES CONNECTIVITES
//==============================================================
void Pyramide::Init_Connectivite()
{
  m_NbVoisins = new int [m_NbSommets];
  m_Voisins   = new int*[m_NbSommets];
  m_NbVoisins[0] = 3;
  m_NbVoisins[1] = 3;
  m_NbVoisins[2] = 3;
  m_NbVoisins[3] = 3;
  m_NbVoisins[4] = 4;
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    m_Voisins  [isommet] = new int[m_NbVoisins[isommet]]; 
 
  m_Voisins[0][0] = 1; m_Voisins[0][1] = 3; m_Voisins[0][2] = 4;
  m_Voisins[1][0] = 2; m_Voisins[1][1] = 0; m_Voisins[1][2] = 4;
  m_Voisins[2][0] = 3; m_Voisins[2][1] = 1; m_Voisins[2][2] = 4;
  m_Voisins[3][0] = 0; m_Voisins[3][1] = 2; m_Voisins[3][2] = 4;
  m_Voisins[4][0] = 0; m_Voisins[4][1] = 1; m_Voisins[4][2] = 2; m_Voisins[4][3] = 3;

  m_SommetsDesAretes = new int*[m_NbAretes];
  for (int iarete = 0; iarete < m_NbAretes; iarete++)
    m_SommetsDesAretes[iarete] = new int[2];

  m_SommetsDesAretes[0][0]=0; m_SommetsDesAretes[0][1]=1;
  m_SommetsDesAretes[1][0]=1; m_SommetsDesAretes[1][1]=2;
  m_SommetsDesAretes[2][0]=2; m_SommetsDesAretes[2][1]=3;
  m_SommetsDesAretes[3][0]=3; m_SommetsDesAretes[3][1]=0;
  m_SommetsDesAretes[4][0]=4; m_SommetsDesAretes[4][1]=0;
  m_SommetsDesAretes[5][0]=4; m_SommetsDesAretes[5][1]=1;
  m_SommetsDesAretes[6][0]=4; m_SommetsDesAretes[6][1]=2;
  m_SommetsDesAretes[7][0]=4; m_SommetsDesAretes[7][1]=3;

  m_NbSommetsDesFaces = new int [m_NbFaces];
  m_SommetsDesFaces   = new int*[m_NbFaces];
  m_NbSommetsDesFaces[0] = 4;
  m_NbSommetsDesFaces[1] = 3;
  m_NbSommetsDesFaces[2] = 3;
  m_NbSommetsDesFaces[3] = 3;
  m_NbSommetsDesFaces[4] = 3;
  for (int i = 0; i < m_NbFaces; i++)
    m_SommetsDesFaces[i]          = new int[m_NbSommetsDesFaces[i]];

  m_SommetsDesFaces[0][0] = 0;  m_SommetsDesFaces[1][0] = 0;  m_SommetsDesFaces[2][0] = 1;  
  m_SommetsDesFaces[0][1] = 1;  m_SommetsDesFaces[1][1] = 1;  m_SommetsDesFaces[2][1] = 2;  
  m_SommetsDesFaces[0][2] = 2;  m_SommetsDesFaces[1][2] = 3;  m_SommetsDesFaces[2][2] = 4;  
  m_SommetsDesFaces[0][3] = 3;                       

  m_SommetsDesFaces[3][0] = 2;  m_SommetsDesFaces[4][0] = 3;    
  m_SommetsDesFaces[3][1] = 3;  m_SommetsDesFaces[4][1] = 0;    
  m_SommetsDesFaces[3][2] = 4;  m_SommetsDesFaces[4][2] = 4;  
}

//==============================================================
//   INITIALISATION DES SOMMETS
//==============================================================
void Pyramide::Init_Sommets(Vecteur *points)
{
  m_Sommets   = new Vecteur[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    m_Sommets[isommet] = points[isommet]; 
}
void Pyramide::Init_SommetsIdeaux()
{
  m_SommetsIdeaux = new Vecteur[m_NbSommets];

  m_SommetsIdeaux[0] =  Vecteur(0,0,0); 
  m_SommetsIdeaux[1] =  Vecteur(1,0,0); 
  m_SommetsIdeaux[2] =  Vecteur(1,1,0); 
  m_SommetsIdeaux[3] =  Vecteur(0,1,0); 
  m_SommetsIdeaux[4] =  Vecteur(.5, .5, .5*std::sqrt(3.)); 

}

//==============================================================
//  JACOBIENNES AUX SOMMETS 
//==============================================================
Matrice Pyramide::Jacobienne(int NumSommet,Vecteur *Noeuds)
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
      J.InitMatrice(2,Noeuds[4]-Noeuds[0]); 
      break;	      
    case 2:	      
      J.InitMatrice(0,Noeuds[2]-Noeuds[3]);
      J.InitMatrice(1,Noeuds[2]-Noeuds[1]); 
      J.InitMatrice(2,Noeuds[4]+Noeuds[2]-Noeuds[1]-Noeuds[3]); 
      break;	      
    case 3:	      
      J.InitMatrice(0,Noeuds[2]-Noeuds[3]);
      J.InitMatrice(1,Noeuds[3]-Noeuds[0]); 
      J.InitMatrice(2,Noeuds[4]-Noeuds[0]); 
      break;
    case 4:	      
      J.InitMatrice(0,Noeuds[2]-Noeuds[3]);
      J.InitMatrice(1,Noeuds[3]-Noeuds[0]); 
      J.InitMatrice(2,Noeuds[4]-Noeuds[0]); 
      break;
    default:
      break;
    }     
  return J;
}

//==============================================================
//  VOLUME
//==============================================================
double Pyramide::VolumeSurface()
{  
  double Resu = 0.;
  int NumSommet;
  int i0,i1,i2;
  Vecteur u, v, w;
  int NbSommetsBase = m_NbSommetsDesFaces[0];

  for (int isommet = 0; isommet < NbSommetsBase; isommet++)
    { 
      NumSommet = m_SommetsDesFaces[0][isommet];
      i0 = m_Voisins[NumSommet][0];
      i1 = m_Voisins[NumSommet][1];
      i2 = m_Voisins[NumSommet][2];
      u  = m_Sommets[i0] - m_Sommets[NumSommet];
      v  = m_Sommets[i1] - m_Sommets[NumSommet];
      w  = m_Sommets[i2] - m_Sommets[NumSommet];
      Resu += u*(v^w)/6;
    }
		      
  return .5*Resu;
}
//==============================================================
//  MESURE D'ODDY MODIFIEE
//==============================================================
double Pyramide::Oddy()
{
  double *D;
  double Dmax;
  Matrice A  (m_Dimension,m_Dimension);
  Matrice Wm1(m_Dimension,m_Dimension);
  Matrice T  (m_Dimension,m_Dimension);
  
  D = new double[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {
      A      = Jacobienne(isommet,m_Sommets);
      Wm1    = Jacobienne(isommet,m_SommetsIdeaux).inverse(); 
      T      = A*Wm1;
      double Det = T.Determinant();
      if (Det == 0)
	{
	  delete [] D;
	  return Infini;
	} 
      else
	{
	  double D1 = std::pow(std::pow(Det,4),1./m_Dimension);
	  double D2 = std::pow((T.transpose()*T).normeFrobenius(),2);
	  double D3 = std::pow(T.normeFrobenius(),4)/m_Dimension;
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
double Pyramide::Conditionnement()
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
//  KNUPP_SKEW : MESURE D'OBLIQUITE
//==============================================================
double Pyramide::Knupp_Skew()
{
  double *kappa;
  double kappaMax;    
  kappa = new double[m_NbSommets];

  Matrice A   (m_Dimension,m_Dimension);

  Matrice Qwm1(m_Dimension,m_Dimension);
  Matrice Q   (m_Dimension,m_Dimension);
  Matrice X   (m_Dimension,m_Dimension);

  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {
      A  = Jacobienne(isommet,m_Sommets);
      if (A.Determinant() == 0)
	{
	  delete [] kappa;
	  return 0.;
	}
      else
	{
	  Qwm1 = Jacobienne(isommet,m_SommetsIdeaux).Decompose3D_Q().inverse(); 
	  Q    = A.Decompose3D_Q();
	  X    = Q*Qwm1;
	  kappa[isommet] = X.normeFrobenius()*(X.inverse().normeFrobenius());
	}
    }
  kappaMax = *PlusGrandElement(kappa,m_NbSommets);
  delete [] kappa;
  return (m_Dimension/ kappaMax);
}

//==============================================================
//   KNUPP_SHAPE : MESURE DE FORME
//==============================================================
double Pyramide::Knupp_Shape()
{
  double *kappa;
  double kappaMax;      
  Matrice A  (m_Dimension,m_Dimension);
  Matrice Wm1(m_Dimension,m_Dimension);
  Matrice T  (m_Dimension,m_Dimension);
  kappa = new double[m_NbSommets];

  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {
      A      = Jacobienne(isommet,m_Sommets);
      Wm1    = Jacobienne(isommet,m_SommetsIdeaux).inverse(); 
      T      = A*Wm1;      
      if (T.Determinant() == 0)
	{
	  delete [] kappa;
	  return 0.;
	}
      else
	kappa[isommet] = T.normeFrobenius()*(T.inverse().normeFrobenius());
    }
  kappaMax = *PlusGrandElement(kappa,m_NbSommets);
  delete [] kappa;
  return (m_Dimension/ kappaMax);
}

//==============================================================
//   KNUPP_VOLUME : MESURE DE VOLUME
//==============================================================
double Pyramide::Knupp_Volume()
{
  double *tau;
  double SommeTau = 0;
  Matrice A  (m_Dimension,m_Dimension);
  Matrice Wm1(m_Dimension,m_Dimension);
  Matrice T  (m_Dimension,m_Dimension);
  
  tau = new double[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {
      A    = Jacobienne(isommet,m_Sommets);
      Wm1  = Jacobienne(isommet,m_SommetsIdeaux).inverse(); 
      T    = A*Wm1;      
      double Det = T.Determinant();
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
double Pyramide::Knupp_VolumeShape()
{
  double *tau;
  double *kappa;
  double kappaMax;
  double SommeTau = 0;
  
  Matrice A  (m_Dimension,m_Dimension);
  Matrice Wm1(m_Dimension,m_Dimension);
  Matrice T  (m_Dimension,m_Dimension);
  
  tau   = new double[m_NbSommets];
  kappa = new double[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {
      A    = Jacobienne(isommet,m_Sommets);
      Wm1  = Jacobienne(isommet,m_SommetsIdeaux).inverse(); 
      T    = A*Wm1;      
      double Det = T.Determinant();
      if ( Det == 0)
	{
	  tau  [isommet] = 0.;
	  kappa[isommet] = Infini;
	}
      else
	{
	  tau[isommet] = (Det>0)?min(Det,1./Det):max(Det,1./Det);
	  kappa[isommet] = T.normeFrobenius()*(T.inverse().normeFrobenius()); 
	} 
      SommeTau +=tau[isommet];      
    }
  kappaMax = *PlusGrandElement(kappa,m_NbSommets);
  delete [] tau;
  delete [] kappa;
  return (SommeTau/m_NbSommets)*(m_Dimension/kappaMax) ;
}

//==============================================================
//  MINIMUM DES JACOBIENS AUX SOMMETS
//==============================================================
double Pyramide::JacobienMin()
{
  double *j;
  double jmin;
  j = new double[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    j[isommet] = Jacobienne(isommet,m_Sommets).Determinant();

  jmin = *PlusPetitElement(j,m_NbSommets);
  delete [] j;
  return jmin;
}

//==============================================================
// GONDOLEMENT DE LA BASE
//==============================================================
double Pyramide::WarpBase()
{
  Vecteur *SommetsBase;
  double Resultat;
  int NbSommetsBase = m_NbSommetsDesFaces[0];
  SommetsBase = new Vecteur[NbSommetsBase];
  for (int i = 0; i < NbSommetsBase; i++)
    SommetsBase[i] = m_Sommets[m_SommetsDesFaces[0][i]];
  Quadrangle Base(m_Dimension, SommetsBase);
  delete [] SommetsBase;
  return Base.Warp();
}
}
