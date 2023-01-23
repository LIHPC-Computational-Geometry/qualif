
#include "Lima/lima++.h"
#include <iostream>

#include "Qualif.h"

//***************************************************************************
//                    EXEMPLE D'UTILISATION DES FONCTIONS
//                            DE LA BIBLIOTHEQUE
//                     DE MESURE DE QUALITE DE MAILLAGE
//***************************************************************************

double const  * min_element(double const * from, double const * end)
{
  double const * min , * current ;
  min = current = from;
  while (++current != end)
    if (*current < * min)
      min = current;
  return min;
}

double const  * max_element(double const * from, double const * end)
{
  double const * max , * current ;
  max = current = from;
  while (++current != end)
    if (*current > * max)
      max = current;
  return max;
}


int main()
{
  Lima::Maillage maillage;
  int NumCritere;
  double *Resultats;
  double const  *ValeurMin, *ValeurMax;
  Qualif::Critere MonCritere;

//==============================================================
//    LECTURE DU MAILLAGE
//==============================================================
  maillage.lire("tests/Maillage3D.unf");

//==============================================================
//    APPEL A AFFICHECRITERES()
//==============================================================
  Qualif::AfficheCriteres();

//==============================================================
//    CHOIX DU CRITERE A APPLIQUER
//==============================================================
  std::cout << " -> NUMERO DU CRITERE A APPLIQUER ? " << std::endl;
  std::cin >> NumCritere;
  std::cout << "CRITERE : " << Qualif::CRITERESTR[NumCritere] <<std::endl; 
  MonCritere = Qualif::Critere(NumCritere);

//==============================================================
//    MESURE DE LA QUALITE DES ELEMENTS
//    DE TOUT LE MAILLAGE : APPELS A ANALYSEPOLYEDRE()
//==============================================================
  std::cout << "" <<std::endl; 
  std::cout << "Mesure de la qualite des polyedres de tout le maillage" <<std::endl;

  int NbElts    = maillage.nb_polyedres();
  Resultats     = new double[NbElts];
  for (int ielement = 0; ielement < NbElts; ielement++)
    Qualif::AnalysePolyedre (maillage.polyedre(ielement),
                               MonCritere,
                               Resultats[ielement]);
  ValeurMin     = min_element(Resultats, Resultats+NbElts);
  ValeurMax     = max_element(Resultats, Resultats+NbElts);
  int NumeroMin = maillage.polyedre(ValeurMin-Resultats).id();
  int NumeroMax = maillage.polyedre(ValeurMax-Resultats).id();
  std::cout << "-> Plus petite valeur du critere : " << *ValeurMin << std::endl;
  std::cout << "   Pour la maille numero         : " <<  NumeroMin << std::endl;
  std::cout << "-> Plus grande valeur du critere : " << *ValeurMax << std::endl;
  std::cout << "   Pour la maille numero         : " <<  NumeroMax <<std::endl;
  delete [] Resultats;
//==============================================================
//    MESURE DE LA QUALITE DES ELEMENTS
//    DES SURFACES DU MAILLAGE : APPELS A ANALYSESURFACE()
//==============================================================
  int NbSurfaces        = maillage.nb_surfaces();
  int DimensionMaillage = maillage.dimension();
    for (int isurface = 0; isurface < NbSurfaces; isurface++)
    {
      int NbEltsSurface    = maillage.surface(isurface).nb_polygones();
      Resultats            = new double[NbEltsSurface];
      std::cout << "" <<std::endl; 
      std::cout << "Mesure  de la qualite des polygones de la surface " 
           << maillage.surface(isurface).nom() <<std::endl; 
      Qualif::AnalyseSurface(maillage.surface(isurface),
                               DimensionMaillage,
                               MonCritere,
                               Resultats);
      ValeurMin = min_element(Resultats, Resultats+NbEltsSurface);
      ValeurMax = max_element(Resultats, Resultats+NbEltsSurface);
      int NumeroMin = maillage.surface(isurface).polygone(ValeurMin-Resultats).id();
      int NumeroMax = maillage.surface(isurface).polygone(ValeurMax-Resultats).id();
      std::cout << "-> Plus petite valeur du critere : " << *ValeurMin << std::endl;
      std::cout << "   Pour la maille numero         : " <<  NumeroMin << std::endl;
      std::cout << "-> Plus grande valeur du critere : " << *ValeurMax << std::endl;
      std::cout << "   Pour la maille numero         : " <<  NumeroMax <<std::endl;
      delete [] Resultats;
    }
}
