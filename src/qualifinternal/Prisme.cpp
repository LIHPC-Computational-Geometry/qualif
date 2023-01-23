#include "Prisme.h"

namespace Qualif {

//**************************************************************
//             METHODES ET DEFINITIONS DE LA CLASSE 
//                            PRISME
//**************************************************************

//==============================================================
//   CONSTRUCTEUR
//==============================================================
Prisme::Prisme(Vecteur *points) 
{  
  m_Dimension      = 3;
  m_NbSommets      = 6;
  m_NbAretes       = 9;
  m_NbFaces        = 5;

  Init_Connectivite();
  Init_Sommets(points);
  Init_SommetsIdeaux();
}
Prisme::Prisme()
{
  m_Dimension      = 3;
  m_NbSommets      = 6;
  m_NbAretes       = 9;
  m_NbFaces        = 5;

  Init_Connectivite();
  Init_SommetsIdeaux();
  m_Sommets   = new Vecteur[m_NbSommets];
}
//==============================================================
//   INITIALISATION DES CONNECTIVITES
//==============================================================
void Prisme::Init_Connectivite()
{
  m_NbVoisins = new int [m_NbSommets];
  m_Voisins   = new int*[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {    
      m_NbVoisins[isommet] = 3;
      m_Voisins  [isommet] = new int[m_NbVoisins[isommet]];
    }
  
  m_Voisins[0][0] = 1; m_Voisins[0][1] = 2; m_Voisins[0][2] = 3;
  m_Voisins[1][0] = 2; m_Voisins[1][1] = 0; m_Voisins[1][2] = 4;
  m_Voisins[2][0] = 0; m_Voisins[2][1] = 1; m_Voisins[2][2] = 5;
  m_Voisins[3][0] = 0; m_Voisins[3][1] = 5; m_Voisins[3][2] = 4;
  m_Voisins[4][0] = 1; m_Voisins[4][1] = 3; m_Voisins[4][2] = 5;
  m_Voisins[5][0] = 2; m_Voisins[5][1] = 4; m_Voisins[5][2] = 3;


  m_SommetsDesAretes = new int*[m_NbAretes];
  for (int iarete = 0; iarete < m_NbAretes; iarete++)
    m_SommetsDesAretes[iarete] = new int[2];

  m_SommetsDesAretes[0][0]=0; m_SommetsDesAretes[0][1]=1;
  m_SommetsDesAretes[1][0]=1; m_SommetsDesAretes[1][1]=2;
  m_SommetsDesAretes[2][0]=2; m_SommetsDesAretes[2][1]=0;
  m_SommetsDesAretes[3][0]=3; m_SommetsDesAretes[3][1]=4;
  m_SommetsDesAretes[4][0]=4; m_SommetsDesAretes[4][1]=5;
  m_SommetsDesAretes[5][0]=5; m_SommetsDesAretes[5][1]=3;
  m_SommetsDesAretes[6][0]=0; m_SommetsDesAretes[6][1]=3;
  m_SommetsDesAretes[7][0]=1; m_SommetsDesAretes[7][1]=4;
  m_SommetsDesAretes[8][0]=2; m_SommetsDesAretes[8][1]=5;

  m_NbSommetsDesFaces= new int [m_NbFaces];
  m_SommetsDesFaces         = new int*[m_NbFaces];
  m_NbSommetsDesFaces[0] = 3;
  m_NbSommetsDesFaces[1] = 3;
  m_NbSommetsDesFaces[2] = 4;
  m_NbSommetsDesFaces[3] = 4;
  m_NbSommetsDesFaces[4] = 4;

  for (int i = 0; i < m_NbFaces; i++)
    m_SommetsDesFaces[i] = new int[m_NbSommetsDesFaces[i]];

  m_SommetsDesFaces[0][0] = 0;  m_SommetsDesFaces[1][0] = 3;    
  m_SommetsDesFaces[0][1] = 1;  m_SommetsDesFaces[1][1] = 5;    
  m_SommetsDesFaces[0][2] = 2;  m_SommetsDesFaces[1][2] = 4;
    
  
  m_SommetsDesFaces[2][0] = 0;  m_SommetsDesFaces[3][0] = 1;  m_SommetsDesFaces[4][0] = 0;    
  m_SommetsDesFaces[2][1] = 3;  m_SommetsDesFaces[3][1] = 2;  m_SommetsDesFaces[4][1] = 2;    
  m_SommetsDesFaces[2][2] = 4;  m_SommetsDesFaces[3][2] = 5;  m_SommetsDesFaces[4][2] = 5;    
  m_SommetsDesFaces[2][3] = 1;  m_SommetsDesFaces[3][3] = 4;  m_SommetsDesFaces[4][3] = 3;    
}

//==============================================================
//   INITIALISATION DES SOMMETS
//==============================================================
void Prisme::Init_Sommets(Vecteur *points)
{
  m_Sommets   = new Vecteur[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    m_Sommets[isommet] = points[isommet]; 
}
void Prisme::Init_SommetsIdeaux()
{
  m_SommetsIdeaux = new Vecteur[m_NbSommets];

  m_SommetsIdeaux[0] =  Vecteur(0,0,0); 
  m_SommetsIdeaux[1] =  Vecteur(1,0,0); 
  m_SommetsIdeaux[2] =  Vecteur(.5,.5*std::sqrt(3.),0) ; 
  m_SommetsIdeaux[3] =  Vecteur(0,0,1); 
  m_SommetsIdeaux[4] =  Vecteur(1,0,1); 
  m_SommetsIdeaux[5] =  Vecteur(.5,.5*std::sqrt(3.),1); 
}

//==============================================================
//  JACOBIENNES AUX SOMMETS 
//==============================================================
Matrice Prisme::Jacobienne(int NumSommet,Vecteur *Noeuds)
{
  Matrice J(m_Dimension,m_Dimension);
  switch(NumSommet)
    {
    case 0:
      J.InitMatrice(0,Noeuds[1]-Noeuds[0]);
      J.InitMatrice(1,Noeuds[2]-Noeuds[0]); 
      J.InitMatrice(2,Noeuds[3]-Noeuds[0]); 
      break;	      
    case 1:	      
      J.InitMatrice(0,Noeuds[1]-Noeuds[0]);
      J.InitMatrice(1,Noeuds[2]-Noeuds[0]); 
      J.InitMatrice(2,Noeuds[4]-Noeuds[1]); 
      break;	      
    case 2:	      
      J.InitMatrice(0,Noeuds[1]-Noeuds[0]);
      J.InitMatrice(1,Noeuds[2]-Noeuds[0]); 
      J.InitMatrice(2,Noeuds[5]-Noeuds[2]); 
      break;	      
    case 3:	      
      J.InitMatrice(0,Noeuds[4]-Noeuds[3]);
      J.InitMatrice(1,Noeuds[5]-Noeuds[3]); 
      J.InitMatrice(2,Noeuds[3]-Noeuds[0]); 
      break;
    case 4:	      
      J.InitMatrice(0,Noeuds[4]-Noeuds[3]);
      J.InitMatrice(1,Noeuds[5]-Noeuds[3]); 
      J.InitMatrice(2,Noeuds[4]-Noeuds[1]); 
      break;
    case 5:	      
      J.InitMatrice(0,Noeuds[4]-Noeuds[3]);
      J.InitMatrice(1,Noeuds[5]-Noeuds[3]); 
      J.InitMatrice(2,Noeuds[5]-Noeuds[2]); 
      break;
    default:
      break;
    }     
  return J;
}

//==============================================================
//  VOLUME
//==============================================================
double Prisme::VolumeSurface()
{  
   double Resu = 0.;
   Vecteur u, v, w;
   u = m_Sommets[1] - m_Sommets[0];
   v = m_Sommets[2] - m_Sommets[0];
   w = m_Sommets[5] - m_Sommets[0];
   Resu += u*(v^w)/6;
   u = m_Sommets[5] - m_Sommets[3];
   v = m_Sommets[4] - m_Sommets[3];
   w = m_Sommets[1] - m_Sommets[3];
   Resu += u*(v^w)/6;
   u = m_Sommets[3] - m_Sommets[1];
   v = m_Sommets[5] - m_Sommets[1];
   w = m_Sommets[0] - m_Sommets[1];
   Resu += u*(v^w)/6;


   u = m_Sommets[5] - m_Sommets[3];
   v = m_Sommets[4] - m_Sommets[3];
   w = m_Sommets[2] - m_Sommets[3];
   Resu += u*(v^w)/6;
   u = m_Sommets[1] - m_Sommets[0];
   v = m_Sommets[2] - m_Sommets[0];
   w = m_Sommets[4] - m_Sommets[0];
   Resu += u*(v^w)/6;
   u = m_Sommets[2] - m_Sommets[0];
   v = m_Sommets[3] - m_Sommets[0];
   w = m_Sommets[4] - m_Sommets[0];
   Resu += u*(v^w)/6;

		      
  return .5*Resu;
}
//==============================================================
//  MESURE D'ODDY MODIFIEE
//==============================================================
double Prisme::Oddy()
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
double Prisme::Conditionnement()
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
double Prisme::ScaledJacobian()
{
  double *scaledJ;
  double scaledJmin;
  Matrice A(m_Dimension,m_Dimension);  
  scaledJ = new double[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {  
      A = Jacobienne(isommet,m_Sommets);
      int i0    = m_Voisins[isommet][0];
      int i1    = m_Voisins[isommet][1];
      int i2    = m_Voisins[isommet][2];
      double l0 =(m_Sommets[isommet] - m_Sommets[i0]).norme2();
      double l1 =(m_Sommets[isommet] - m_Sommets[i1]).norme2();
      double l2 =(m_Sommets[isommet] - m_Sommets[i2]).norme2();
      double l012 = l0*l1*l2;
      double Det = A.Determinant();
      if (l012 == 0 || Det == 0)
	scaledJ[isommet] = 0.;
      else
	scaledJ[isommet] = Det/l012;
    }
  scaledJmin = *PlusPetitElement(scaledJ,m_NbSommets);
  delete [] scaledJ;

  scaledJmin *= 2./sqrt(3.);
  if(scaledJmin >  1.) scaledJmin =  1.;
  if(scaledJmin < -1.) scaledJmin = -1.;

  return scaledJmin;
}

//==============================================================
//  KNUPP_SKEW : MESURE D'OBLIQUITE
//==============================================================
double Prisme::Knupp_Skew()
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
      A    = Jacobienne(isommet,m_Sommets);
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
double Prisme::Knupp_Shape()
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
double Prisme::Knupp_Volume()
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
double Prisme::Knupp_VolumeShape()
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
double Prisme::JacobienMin()
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

}


