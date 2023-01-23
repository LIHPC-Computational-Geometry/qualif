#ifndef DEF_QUALMAIL_H
#define DEF_QUALMAIL_H

#include <Lima/lima++.h>
#include "ConfigQualif.h"
#include "ConstQualif.h"
//***************************************************************************
//
//                        FICHIER INTERFACE
//                           DE LA BIBLIOTHEQUE 
//                     DE MESURE DE QUALITE DE MAILLAGE
//
//***************************************************************************

namespace Qualif {

//==============================================================
//   PROTOTYPES DES FONCTIONS UTILISATEURS DE LA BIBLIOTHEQUE
//==============================================================
void AfficheCriteres ();
void AnalyseVolume   (const Lima::Volume   volume,   Critere critere, double *Resultats);
void AnalysePolyedre (const Lima::Polyedre polyedre, Critere critere, double &Resultat);
void AnalyseSurface  (const Lima::Surface  surface, int DimensionMailage, 
		      Critere critere, double *Resultats);
void AnalysePolygone (const Lima::Polygone polygone, int DimensionMailage, 
		      Critere critere, double &Resultat);
const std::string QualifVersion();
}
#endif
