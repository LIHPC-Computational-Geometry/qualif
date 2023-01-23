#include <fstream>


#include "Qualif.h"
#include "Fonctions.h"
#include "EcrireLima.h"
using namespace Qualif;

//***********************************************
// MESURE DE LA QUALITE DES MAILLAGES 
// POUR LA VALIDATION
//

// AFFICHE:
//    LE NOM DU CRITERE
//    LE NB D'ELTS DU MAILLAGE
//    LE NB D'ELTS MESURES
//    LA VALEUR MINIMUM DU CRITERE
//    LA VALEUR MAXIMUM DU CRITERE
//    LA MOYENNE
//    L'ECART TYPE

// CREE:
//   LES MAILLAGES 'NOM_DU_CRITERE'_Max.unf
//              ET 'NOM_DU_CRITERE'_Min.unf
//   RESPECTIVEMENT :
//   LA MAILLE CORRESPONDANT A LA VALEUR MAXIMALE
//   LA MAILLE CORRESPONDANT A LA VALEUR MINIMALE
//
//***********************************************

void statistique (double * Resultats, int NbElts)
{
      int NbEltsMesures = NbElts;

     double Moyenne = 0.0;
     
      for (int ielement = 0; ielement < NbElts; ++ ielement)  
      if (Resultats[ielement] == Resultats[ielement])
        Moyenne += Resultats[ielement];
      else
       NbEltsMesures --;
   Moyenne    = Moyenne/(double)NbEltsMesures;
   double  * ValeurMin  = Qualif::PlusPetitElement(Resultats,NbElts);
   double  * ValeurMax  = Qualif::PlusGrandElement(Resultats,NbElts);
   double EcartType = 0.0;
   
   for (int ielement = 0; ielement < NbElts; ielement++)
     if (Resultats[ielement] == Resultats[ielement])
	  EcartType += std::pow(Resultats[ielement]-Moyenne,2);
   
    EcartType = std::sqrt(EcartType/(double)NbEltsMesures);

      std::cout << " "<<std::endl;
      std::cout << "========================"<<std::endl;
      std::cout <<"Nb Elts du maillage : "<<std::setw(10)<< NbElts<<std::endl;
      std::cout <<"Nb Elts mesures     : "<<std::setw(10)<< NbEltsMesures<<std::endl;
      
   int IndiceMin, IndiceMax;
   double ElementMin, ElementMax;
   
      if (ValeurMin != 0 || ValeurMax != 0) 
	{
	  ElementMin = ValeurMin - Resultats;
	  ElementMax = ValeurMax - Resultats;

	  std::cout <<"Valeur minimum      : "<<std::setw(10)<< *ValeurMin
	       <<" (Element "<<ElementMin<<" )"<<std::endl; 
	  std::cout <<"Valeur maximum      : "<<std::setw(10)<< *ValeurMax
	       <<" (Element "<<ElementMax<<" )"<<std::endl; 
	  std::cout <<"Valeur Moyenne      : "<<std::setw(10)<< Moyenne << std::endl;
	  std::cout <<"Ecart Type          : "<<std::setw(10)<< EcartType << std::endl;    
	}
}

	  
int main(int argc, char * argv[])
{

  char NomMaillage[60];
  Lima::Maillage maya;
  int DimensionMaillage;
  double *Resultats;

  int    NumMaillage;
  int    NumCritere,NumVolume,NumSurface;
  int    NbElts, NbEltsMesures;
  double *ValeurMin, *ValeurMax;
  int    IndiceMin, IndiceMax;
  int    ElementMin, ElementMax;
  double Moyenne   = 0.;
  double EcartType = 0.;
  bool   TraitePolyedre = false;
  char   rep;
  
  std::cout << "     ************************************" << std::endl;
  std::cout << "         MESURE DE QUALITE DE MAILLAGE   " << std::endl;
  std::cout << "     ************************************" << std::endl;
  std::cout << " "<< std::endl;

 
  if (argc != 2)
   {
   	std::cerr << "Usage " 
   	          << argv[0] << " nom_maillage" << std::endl;
	          exit (-1);
   }
   
  //======================================
  //         LECTURE DU MAILLAGE
  //======================================
  maya.lire(argv[1]);  
  DimensionMaillage = maya.dimension();
  if (DimensionMaillage == 3 && maya.nb_polyedres() != 0)
    TraitePolyedre = true;      

  char continuer = 'o';
  while (continuer == 'o' || continuer == 'O')
    {
      Qualif::AfficheCriteres();
      std::cout << " -> NUMERO DU CRITERE A APPLIQUER ? " << std::endl;
      std::cin >> NumCritere;
      std::cout << " "<< std::endl;     
      std::cout << " APPLICATION DU CRITERE : " << Qualif::CRITERESTR[NumCritere] << std::endl;
      
      //======================================
      //        APPEL A Qualif
      //======================================
      if (TraitePolyedre)
	{
	      for (int ivolume = 0; ivolume < maya.nb_volumes(); ivolume++)
	      {
		std::cout << "Volume : " << ivolume << "." << maya.volume(ivolume).nom()<< std::endl;
                NbElts          = maya.volume(ivolume).nb_polyedres();
	        NbEltsMesures   = NbElts;
	        Resultats = new double[NbElts];
	        Qualif::AnalyseVolume (maya.volume(ivolume), 
				       Qualif::Critere(NumCritere),
				       Resultats);

               statistique (Resultats, NbElts);
               delete [] Resultats;

	      }


        }
      else
	{
	  if (maya.nb_surfaces() > 0) 
	    {
	      std::cout << "  Traitement par surface ? o/n" << std::endl;
	      std::cin >> rep;	  
	    }
	  if (rep == 'o' || rep == 'O')
	    {
	      for (int isurface = 0; isurface < maya.nb_surfaces(); isurface++)
		std::cout << isurface << "." << maya.surface(isurface).nom()<< std::endl;
	      std::cout << "Numero de la surface a traiter ?"<< std::endl;
	      std::cin >> NumSurface;
	      NbElts        = maya.surface(NumSurface).nb_polygones();
	      NbEltsMesures = NbElts;
	      Resultats = new double[NbElts];
	      Qualif::AnalyseSurface (maya.surface(NumSurface),
					DimensionMaillage,
					Qualif::Critere(NumCritere),
					Resultats);
	    }
	  else
	    {
	      NbElts        = maya.nb_polygones();
	      NbEltsMesures = NbElts;
	      Resultats     = new double[NbElts];
	      for (int ielement = 0; ielement < NbElts; ielement++)
		{
		  Qualif::AnalysePolygone(maya.polygone(ielement),
					    DimensionMaillage,
					    Qualif::Critere(NumCritere),
					    Resultats[ielement]);

		}	     
	      for (int ielement = 0; ielement < NbElts; ielement++)
		if (Resultats[ielement] == Resultats[ielement])
		  Moyenne += Resultats[ielement];
		else
		  NbEltsMesures --;
		  	
	    }
     
           statistique (Resultats, NbElts);
           delete [] Resultats;
     
	}

  
      std::cout << "Autre Critere ? o/n" << std::endl;
      std::cin >> continuer;  
    }
}



