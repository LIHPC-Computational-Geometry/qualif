#include "IQualif.h"
#include "Qualif.h"
#include "Maille.h"
#include "Triangle.h"
#include "Quadrangle.h"
#include "Tetraedre.h"
#include "Pyramide.h"
#include "Prisme.h"
#include "Hexaedre.h"

//***************************************************************************
//
//                       FONCTIONS UTILISATEUR 
//                         DE LA BIBLIOTHEQUE 
//                             DE MESURE 
//                      DE QUALITE DE MAILLAGE
//
//***************************************************************************
// NB il s'agit de l'interface entre Lima et Qualif(internal) [EB]

namespace Qualif {

//==============================================================
//    AFFICHAGE DES CRITERES DISPONIBLES DANS LA BIBLIOTHEQUE
//==============================================================
void AfficheCriteres()
{
  std::cout << " " << std::endl;
  std::cout << " CRITERES DISPONIBLES DANS LA BIBLIOTHEQUE DE MESURE : " << std::endl; 
  std::cout << " ----------------------------------------------------- " << std::endl;
  std::cout << " " << std::endl;
  for (int icritere = 0; icritere < FIN; icritere++)
    std::cout << " " << std::setw(2) << icritere << ". " << CRITERESTR[icritere] << std::endl;
  std::cout << " " << std::endl;
}

//==============================================================
//    ANALYSE D'UNE SURFACE LIMA PAR L'APPLICATION D'UN CRITERE
//            BOUCLE SUR LES POLYGONES DE LA SURFACE
//==============================================================
void AnalyseSurface(const Lima::Surface surface, int DimensionMaillage, 
		    Critere critere, double *Resultats)
{
  for (int ipolygone = 0; ipolygone < surface.nb_polygones(); ipolygone++)
    AnalysePolygone(surface.polygone(ipolygone),DimensionMaillage,
		    critere,Resultats[ipolygone]);
}

//==============================================================
//    ANALYSE D'UN VOLUME LIMA PAR L'APPLICATION D'UN CRITERE
//            BOUCLE SUR LES POLYEDRES DU VOLUME
//==============================================================
void AnalyseVolume(const Lima::Volume volume, Critere critere, double *Resultats)
{
  for (int ipolyedre = 0; ipolyedre < volume.nb_polyedres(); ipolyedre++)
    AnalysePolyedre(volume.polyedre(ipolyedre),critere,Resultats[ipolyedre]);
}

//==============================================================
//    ANALYSE D'UN POLYGONE PAR L'APPLICATION D'UN CRITERE
//==============================================================
void AnalysePolygone(const Lima::Polygone polygone, int DimensionMaillage,
		      Critere critere, double &Resultat)
{
  int nb_sommets;
  Vecteur *sommets;
  Maille  *MailleCourante;

  nb_sommets = polygone.nb_noeuds();
  sommets    = new Vecteur[nb_sommets];

  for (int i=0;i<nb_sommets;i++)
    sommets[i] = Vecteur(polygone.noeud(i).x(),
			 polygone.noeud(i).y(),
			 polygone.noeud(i).z());

  if (nb_sommets == 3)
    {
      MailleCourante = new Triangle(DimensionMaillage,sommets);
      Resultat = MailleCourante->AppliqueCritere(critere);
      //MailleCourante->Detruit();
    }
   else if (nb_sommets == 4)
    {
      MailleCourante = new Quadrangle(DimensionMaillage,sommets); 
      Resultat = MailleCourante->AppliqueCritere(critere);
      //MailleCourante->Detruit();
    }
  else
    {
      std::cout << "*** NOMBRE DE SOMMETS NON PRIS EN COMPTE ***" << std::endl;
      Resultat = Not_A_Number;
    }      
  delete MailleCourante;
  delete[] sommets;
}

//==============================================================
//    ANALYSE D'UN POLYEDRE PAR L'APPLICATION D'UN CRITERE
//==============================================================
void AnalysePolyedre(const Lima::Polyedre polyedre, Critere critere, double &Resultat)
{
  int nb_sommets;
  Vecteur *sommets;
  Maille  *MailleCourante;

  nb_sommets = polyedre.nb_noeuds();
  sommets    = new Vecteur[nb_sommets];

  for (int i=0;i<nb_sommets;i++)
    sommets[i] = Vecteur(polyedre.noeud(i).x(),
			 polyedre.noeud(i).y(),
			 polyedre.noeud(i).z());
  
  switch (polyedre.type()) 
    {
    case Lima::Polyedre::TETRAEDRE : 
      MailleCourante = new Tetraedre(sommets);
      Resultat       = MailleCourante->AppliqueCritere(critere);
      //MailleCourante->Detruit();
      break;
    case Lima::Polyedre::PYRAMIDE :
      MailleCourante = new Pyramide(sommets);
      Resultat       = MailleCourante->AppliqueCritere(critere);
      //MailleCourante->Detruit();
      break;
    case Lima::Polyedre::PRISME :
      MailleCourante = new Prisme(sommets);
      Resultat       = MailleCourante->AppliqueCritere(critere);
      //MailleCourante->Detruit();
      break;
    case Lima::Polyedre::HEXAEDRE :
      MailleCourante =  new Hexaedre(sommets);
      Resultat       = MailleCourante->AppliqueCritere(critere);
      //MailleCourante->Detruit();
      break;
    default :
      std::cout << "***    TYPE NON IMPLEMENTE...     ***" << std::endl;
      Resultat = Not_A_Number;
      break;
    }
  delete MailleCourante;
  delete[] sommets;
}
  
}
