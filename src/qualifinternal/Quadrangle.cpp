#include "Quadrangle.h"

namespace Qualif {

//**************************************************************
//             METHODES ET DEFINITIONS DE LA CLASSE 
//                            QUADRANGLE
//**************************************************************
//==============================================================
//   CONSTRUCTEUR
//==============================================================
Quadrangle::Quadrangle(int DimMail, Vecteur *points) 
{  
  m_DimensionMaillage = DimMail;

  m_Dimension      = 2;
  m_NbSommets      = 4;
  m_NbAretes       = 4;

  Init_Connectivite();
  Init_Sommets(points);
}
Quadrangle::Quadrangle(int DimMail)
{
  m_DimensionMaillage = DimMail;

  m_Dimension      = 2;
  m_NbSommets      = 4;
  m_NbAretes       = 4;

  Init_Connectivite();
  m_Sommets   = new Vecteur[m_NbSommets];
}
//==============================================================
//   INITIALISATION DES CONNECTIVITES
//==============================================================
void Quadrangle::Init_Connectivite()
{
  m_NbVoisins = new int [m_NbSommets];
  m_Voisins   = new int*[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {    
      m_NbVoisins[isommet] = 2;
      m_Voisins  [isommet] = new int[m_NbVoisins[isommet]];
    }
  
  m_Voisins[0][0] = 1; m_Voisins[0][1] = 3;
  m_Voisins[1][0] = 2; m_Voisins[1][1] = 0;
  m_Voisins[2][0] = 3; m_Voisins[2][1] = 1;
  m_Voisins[3][0] = 0; m_Voisins[3][1] = 2;

  m_SommetsDesAretes = new int*[m_NbAretes];
  for (int iarete = 0; iarete < m_NbAretes; iarete++)
    m_SommetsDesAretes[iarete] = new int[2];

  m_SommetsDesAretes[0][0]=0; m_SommetsDesAretes[0][1]=1;
  m_SommetsDesAretes[1][0]=1; m_SommetsDesAretes[1][1]=2;
  m_SommetsDesAretes[2][0]=2; m_SommetsDesAretes[2][1]=3;
  m_SommetsDesAretes[3][0]=3; m_SommetsDesAretes[3][1]=0;
}

//==============================================================
//   INITIALISATION DES SOMMETS
//==============================================================
void Quadrangle::Init_Sommets(Vecteur *points)
{
  m_Sommets   = new Vecteur[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    m_Sommets[isommet] = points[isommet]; 
}

//==============================================================
//  JACOBIENNES AUX SOMMETS POUR UN QUADRANGLE
//==============================================================
Matrice Quadrangle::Jacobienne(int NumSommet,Vecteur *Noeuds)
{
  Matrice J(m_Dimension,m_Dimension);

  if (m_DimensionMaillage==2)
    switch(NumSommet)
      {
      case 0:
	J.InitMatrice(0,Noeuds[1]-Noeuds[0]);
	J.InitMatrice(1,Noeuds[3]-Noeuds[0]); 
	break;	      	      		
      case 1:	      	      		
	J.InitMatrice(0,Noeuds[1]-Noeuds[0]);
	J.InitMatrice(1,Noeuds[2]-Noeuds[1]); 
	break;	      	      		
      case 2:	      	      		
	J.InitMatrice(0,Noeuds[2]-Noeuds[3]);
	J.InitMatrice(1,Noeuds[2]-Noeuds[1]); 
	break;	      	      		
      case 3:	      	      		
	J.InitMatrice(0,Noeuds[2]-Noeuds[3]);
	J.InitMatrice(1,Noeuds[3]-Noeuds[0]); 
	break;
      }  
  else if (m_DimensionMaillage==3)
    {
      Vecteur *NoeudsNouveauRepere = new Vecteur[3];
      int Num0 = m_Voisins[NumSommet][0];
      int Num1 = m_Voisins[NumSommet][1];
      Noeuds[NumSommet].NouveauRepere(Noeuds[Num0],Noeuds[Num1],
				      NoeudsNouveauRepere);
      J.InitMatrice(0,NoeudsNouveauRepere[1]-NoeudsNouveauRepere[0]);
      J.InitMatrice(1,NoeudsNouveauRepere[2]-NoeudsNouveauRepere[0]);   
      delete [] NoeudsNouveauRepere;
    }
  return J;
}
//==============================================================
//  JACOBIENNE AU CENTRE
//==============================================================
Matrice Quadrangle::JacobienneAuCentre(Vecteur *Noeuds)
{

  Vecteur u[2];
  if (m_DimensionMaillage == 2)
    {
      Vecteur Nul(0,0);
      u[0] = (Nul-Noeuds[0]+Noeuds[1]+Noeuds[2]-Noeuds[3])*.5;
      u[1] = (Nul-Noeuds[0]-Noeuds[1]+Noeuds[2]+Noeuds[3])*.5;
    }
  else if (m_DimensionMaillage == 3)
    {
      Vecteur Centre = (m_Sommets[0]+m_Sommets[1]+
			m_Sommets[2]+m_Sommets[3])*0.25;
      Vecteur MilieuArete0 = (m_Sommets[0]+m_Sommets[1])*.5;
      Vecteur MilieuArete1 = (m_Sommets[1]+m_Sommets[2])*.5;

      Vecteur *NoeudsNouveauRepere = new Vecteur[3];
      Centre.NouveauRepere(MilieuArete0,MilieuArete1,
			   NoeudsNouveauRepere);
      u[0] = (NoeudsNouveauRepere[1] - NoeudsNouveauRepere[0])*2.;
      u[1] = (NoeudsNouveauRepere[2] - NoeudsNouveauRepere[0])*2.;
      delete [] NoeudsNouveauRepere;
   }      
      return Matrice(m_Dimension,m_Dimension,u);
}
//==============================================================
//  SURFACE
//==============================================================
double Quadrangle::VolumeSurface()
{  
  double Resu = 0;
  Vecteur u, v;
  Vecteur Centre = (m_Sommets[0]+m_Sommets[1]+
		    m_Sommets[2]+m_Sommets[3])*0.25;
  
  for (int iarete = 0; iarete < m_NbAretes; iarete++)
    {      
      int i0 = m_SommetsDesAretes[iarete][0];
      int i1 = m_SommetsDesAretes[iarete][1];
      if (m_DimensionMaillage == 2)
	{
	  u = m_Sommets[i0] - Centre;
	  v = m_Sommets[i1] - Centre;
	  Resu += 0.5*u.Determinant(v);
	}
      else if (m_DimensionMaillage == 3)
	{
	  Vecteur *NoeudsNouveauRepere = new Vecteur[3];
	  Centre.NouveauRepere(m_Sommets[i0],m_Sommets[i1],
				   NoeudsNouveauRepere);
	  u = NoeudsNouveauRepere[1] - NoeudsNouveauRepere[0];
	  v = NoeudsNouveauRepere[2] - NoeudsNouveauRepere[0];
	  delete [] NoeudsNouveauRepere;
//	  std::cout << 0.5*u.Determinant(v)<< std::endl;
	  Resu += 0.5*u.Determinant(v);
	}      
    }
  return Resu;
}

//==============================================================
//  MESURE D'ODDY MODIFIEE
//==============================================================
double Quadrangle::Oddy()
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
double Quadrangle::Conditionnement()
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
double Quadrangle::ScaledJacobian()
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
      double l0 =(m_Sommets[isommet] - m_Sommets[i0]).norme2();
      double l1 =(m_Sommets[isommet] - m_Sommets[i1]).norme2();
      double l01 = l0*l1;
      double Det = A.Determinant();
      if (l01 == 0 || Det == 0)
	scaledJ[isommet] = 0.;
      else
	scaledJ[isommet] = Det/l01;
    }
  scaledJmin = *PlusPetitElement(scaledJ,m_NbSommets);
  delete [] scaledJ;
  return scaledJmin;
}

//==============================================================
//  KNUPP_SKEW : MESURE D'OBLIQUITE
//==============================================================
double Quadrangle::Knupp_Skew()
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
	  Q = A.Decompose2D_Q();
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
double Quadrangle::Knupp_Shape()
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
double Quadrangle::Knupp_Volume()
{
  double *tau;
  double SommeTau = 0;
  Matrice A(m_Dimension,m_Dimension);
  
  tau = new double[m_NbSommets];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    {
      A = Jacobienne(isommet,m_Sommets);
      double Det = A.Determinant();
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
double Quadrangle::Knupp_VolumeShape()
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
//   ANGLE MINIMUM OU MAXIMUM
//==============================================================
double Quadrangle::AngleMinMax(const std::string MinMax)
{

  double *angle;
  double angleminmax;  
  Vecteur u,v;

  if (VolumeSurface()==0) return 0.;

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
  return 180./PI* (angleminmax); 
}

//==============================================================
//  MINIMUM DES JACOBIENS AUX SOMMETS ET AU CENTRE 
//==============================================================
double Quadrangle::JacobienMin()
{
  double *j;
  double jmin;
  j = new double[m_NbSommets+1];
  for (int isommet = 0; isommet < m_NbSommets; isommet++)
    j[isommet] = Jacobienne(isommet,m_Sommets).Determinant();

  j[m_NbSommets] = JacobienneAuCentre(m_Sommets).Determinant();

  jmin = *PlusPetitElement(j,m_NbSommets+1);
  delete [] j;
  return jmin;
}

//==============================================================
//  JACOBIEN AU CENTRE 
//==============================================================
double Quadrangle::JacobienCentre()
{
  return JacobienneAuCentre(m_Sommets).Determinant();
}

//==============================================================
//   RAPPORT D'ASPECT
//==============================================================
double Quadrangle::AspectRatioCenter()
{
  Vecteur MilieuArete1, MilieuArete2;
  Vecteur Centre;
  Vecteur a, b;
      
  MilieuArete1 = (m_Sommets[1] + m_Sommets[2])/2.;
  MilieuArete2 = (m_Sommets[2] + m_Sommets[3])/2.;
  Centre       = (m_Sommets[0] + m_Sommets[1] + 
		  m_Sommets[2] + m_Sommets[3])/4.;
      
  a  = MilieuArete1-Centre;
  b  = MilieuArete2-Centre;
      
  double l1 = a.norme2();
  double l2 = b.norme2();

  if (l1>l2)
    return l1/l2;
  else
    return l2/l1;
}

//==============================================================
//   MESURE D'OBLIQUITE 
//==============================================================
double Quadrangle::Skew()
{
  Vecteur MilieuArete1, MilieuArete2;
  Vecteur Centre;
  Vecteur a1,a2;
  double  l1,l2;
  
  MilieuArete1 = (m_Sommets[1] + m_Sommets[2])/2.;
  MilieuArete2 = (m_Sommets[2] + m_Sommets[3])/2.;
  Centre       = (m_Sommets[0] + m_Sommets[1] + 
		  m_Sommets[2] + m_Sommets[3])/4.;
  
  a1  = MilieuArete1-Centre;
  a2  = MilieuArete2-Centre;
  l1  = a1.norme2();
  l2  = a2.norme2();

  if (l1 == 0 ||l2 == 0)
    return 1.;
  else
    return std::fabs((a1*a2)/(l1*l2));
  
}

//==============================================================
//   TAPERS DE ROBINSON (REVUS PAR FIELD)
//==============================================================
double Quadrangle::Taper()
{
  if (m_DimensionMaillage == 2)
    {
      Vecteur MilieuArete1, MilieuArete2;
      Vecteur MilieuDiag1, MilieuDiag2;
      Vecteur Centre, C;
      Vecteur a, b;
      double psX,psY;
      double tauX,tauY;

      MilieuArete1 = (m_Sommets[1] + m_Sommets[2])*0.5;
      MilieuArete2 = (m_Sommets[2] + m_Sommets[3])*0.5;
      Centre       = (m_Sommets[0] + m_Sommets[1] + 
		      m_Sommets[2] + m_Sommets[3])*0.25;
       
      MilieuDiag1 = (m_Sommets[0] + m_Sommets[2])*0.5;
      MilieuDiag2 = (m_Sommets[1] + m_Sommets[3])*0.5;
      C           = (MilieuDiag1-MilieuDiag2)*0.5;

      a  = MilieuArete1-Centre;
      b  = MilieuArete2-Centre; 
  
      Vecteur a_ortho(-a.y(),a.x());
      psX = b*a_ortho;
      psY = a*a;
      
      if (psX == 0 || psY == 0)
	return 0.;
      else
	{
	  tauX = std::fabs((C*a_ortho)/psX);     
	  tauY = std::fabs((C*a)/psY); 
	  return max(tauX,tauY);
	}
    }
  else
    return Not_A_Number;
}

//==============================================================
//  ETIREMENT 
//==============================================================
double Quadrangle::Etirement()
{
  double *L;
  double Lmin;
  double Diag[2];
  double Diagmax;
  L = new double [m_NbAretes];
  for (int iarete = 0; iarete < m_NbAretes; iarete++)
    {
      int i0 = m_SommetsDesAretes[iarete][0];
      int i1 = m_SommetsDesAretes[iarete][1];
      L[iarete]  = (m_Sommets[i1] - m_Sommets[i0]).norme2();
    }
  Lmin = *PlusPetitElement(L,m_NbAretes);

  Diag[0] = (m_Sommets[2] - m_Sommets[0]).norme2();
  Diag[1] = (m_Sommets[3] - m_Sommets[1]).norme2();
  Diagmax = max(Diag[0],Diag[1]);
  delete [] L;
  return std::sqrt(2.)*Lmin/Diagmax;
}

//==============================================================
//  DEVIATION PAR RAPPORT AU PLAN 
//==============================================================
double Quadrangle::Warp()
{
  if (m_DimensionMaillage == 2)
    return 0.;
  else
    {
      double *L;
      Vecteur *MilieuArete;
      Vecteur Centre;
      Vecteur u1, u2;
      Vecteur u_ortho;
      Vecteur v;
      int IndiceAreteMin;
      double N_u_ortho, N_v;
      L           = new double  [m_NbAretes];
      MilieuArete = new Vecteur [m_NbAretes];
      
      for (int iarete = 0; iarete < m_NbAretes; iarete++)
	{
	  int i0 = m_SommetsDesAretes[iarete][0];
	  int i1 = m_SommetsDesAretes[iarete][1];
	  L[iarete]           = (m_Sommets[i1] - m_Sommets[i0]).norme2();
	  MilieuArete[iarete] = (m_Sommets[i0] + m_Sommets[i1])/2.;
	}
      IndiceAreteMin = PlusPetitElement(L,m_NbAretes) - L;
      
      Centre         = (m_Sommets[0] + m_Sommets[1] + 
		        m_Sommets[2] + m_Sommets[3])*.25;
  
      u1  = MilieuArete[1] - Centre;
      u2  = MilieuArete[2] - Centre;
      
      u_ortho = u1^u2;
      v       = (MilieuArete[IndiceAreteMin] - 
	        m_Sommets[m_SommetsDesAretes[IndiceAreteMin][1]])*.5;
      delete [] L;
      delete [] MilieuArete;
  
      N_u_ortho = u_ortho.norme2();
      N_v       = v.norme2();

      if ( N_u_ortho == 0 || N_v == 0)
	return Infini;
      else
	return 180./PI * std::fabs(std::asin((u_ortho*v)/(N_u_ortho*N_v)));
    } 
}
}

