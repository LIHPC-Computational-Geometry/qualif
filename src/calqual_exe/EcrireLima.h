#include "Qualif.h"
#include <strings.h>



void EcrireLimaPolygone(Lima::Maillage maya, int ElementMin, int ElementMax, int NumCritere)
{
  int    IndiceMin;
  int    IndiceMax;
  int    ElementMinMax;
  char   fichier[40];
  Lima::Maillage mayaMinMax[2];

  for (int iminmax = 0; iminmax < 2; iminmax++)
    {
      switch (iminmax)
	{
	case 0:
	  ElementMinMax = ElementMin; 
	  std::strcpy(fichier,Qualif::CRITERESTR[NumCritere].c_str());
	  std::strcat(fichier,"_Min.unf");
	  break;
	case 1:
	  ElementMinMax = ElementMax; 
	  std::strcpy(fichier,Qualif::CRITERESTR[NumCritere].c_str());
	  std::strcat(fichier,"_Max.unf");
	  break;
	default:
	  break;
	}   
      mayaMinMax[iminmax].dimension( maya.dimension());
      for (int isommet = 0; isommet < maya.polygone(ElementMinMax).nb_noeuds(); isommet++)
	{
	  Lima::Noeud n(maya.polygone(ElementMinMax).noeud(isommet).id(),
			maya.polygone(ElementMinMax).noeud(isommet).x(),
			maya.polygone(ElementMinMax).noeud(isommet).y(),
			maya.polygone(ElementMinMax).noeud(isommet).z());
	  mayaMinMax[iminmax].ajouter(n);
	}
      switch (maya.polygone(ElementMinMax).nb_noeuds())
	{
	case 3:
	  {
	    Lima::Noeud n0 = mayaMinMax[iminmax].noeud_id(maya.polygone(ElementMinMax).noeud(0).id()); 
	    Lima::Noeud n1 = mayaMinMax[iminmax].noeud_id(maya.polygone(ElementMinMax).noeud(1).id()); 
	    Lima::Noeud n2 = mayaMinMax[iminmax].noeud_id(maya.polygone(ElementMinMax).noeud(2).id()); 
	    Lima::Polygone poly (maya.polygone(ElementMinMax).id(),
				 n0,n1,n2);
	    mayaMinMax[iminmax].ajouter(poly);
	  }
	  break;
	case 4:
	  {
	    Lima::Noeud n0 = mayaMinMax[iminmax].noeud_id(maya.polygone(ElementMinMax).noeud(0).id()); 
	    Lima::Noeud n1 = mayaMinMax[iminmax].noeud_id(maya.polygone(ElementMinMax).noeud(1).id()); 
	    Lima::Noeud n2 = mayaMinMax[iminmax].noeud_id(maya.polygone(ElementMinMax).noeud(2).id()); 
	    Lima::Noeud n3 = mayaMinMax[iminmax].noeud_id(maya.polygone(ElementMinMax).noeud(3).id()); 
	    Lima::Polygone poly (maya.polygone(ElementMinMax).id(),
				 n0,n1,n2,n3);
	    mayaMinMax[iminmax].ajouter(poly);
	  }
	  break;
	default:
	  break;
	}
      mayaMinMax[iminmax].ecrire(fichier);
    }
}


void EcrireLimaPolyedre(Lima::Maillage maya, int ElementMin, int ElementMax, int NumCritere)
{
  int    IndiceMin;
  int    IndiceMax;
  int    ElementMinMax;
  char   fichier[40];
  Lima::Maillage mayaMinMax[2];

  for (int iminmax = 0; iminmax < 2; iminmax++)
    {
      switch (iminmax)
	{
	case 0:
	  ElementMinMax = ElementMin; 
	  std::strcpy(fichier,Qualif::CRITERESTR[NumCritere].c_str());
	  std::strcat(fichier,"_Min.unf");
	  break;
	case 1:
	  ElementMinMax = ElementMax; 
	  std::strcpy(fichier,Qualif::CRITERESTR[NumCritere].c_str());
	  std::strcat(fichier,"_Max.unf");
	  break;
	default:
	  break;
	}   
      
      mayaMinMax[iminmax].dimension(maya.dimension());
      Lima::Volume vol("VOLUME");
      for (int isommet = 0; isommet < maya.polyedre(ElementMinMax).nb_noeuds(); isommet++)
	{
	  Lima::Noeud n(maya.polyedre(ElementMinMax).noeud(isommet).id(),
			maya.polyedre(ElementMinMax).noeud(isommet).x(),
			maya.polyedre(ElementMinMax).noeud(isommet).y(),
			maya.polyedre(ElementMinMax).noeud(isommet).z());
	  mayaMinMax[iminmax].ajouter(n);
	}
      switch (maya.polyedre(ElementMinMax).nb_noeuds())
	{
	case 4:
	  {
	    Lima::Noeud n0 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(0).id()); 
	    Lima::Noeud n1 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(1).id()); 
	    Lima::Noeud n2 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(2).id()); 
	    Lima::Noeud n3 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(3).id()); 
	    Lima::Polyedre poly (maya.polyedre(ElementMinMax).id(),n0,n1,n2,n3);
	    mayaMinMax[iminmax].ajouter(poly);
	    mayaMinMax[iminmax].ajouter(vol);
	    vol.ajouter(poly);
	  }
	  break;
	case 5:
	  {
	    Lima::Noeud n0 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(0).id()); 
	    Lima::Noeud n1 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(1).id()); 
	    Lima::Noeud n2 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(2).id()); 
	    Lima::Noeud n3 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(3).id()); 
	    Lima::Noeud n4 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(4).id()); 
	    Lima::Polyedre poly (maya.polyedre(ElementMinMax).id(),n0,n1,n2,n3,n4);
	    mayaMinMax[iminmax].ajouter(poly);
	    mayaMinMax[iminmax].ajouter(vol);
	    vol.ajouter(poly);
	  }
	  break;
	case 6:
	  {
	    Lima::Noeud n0 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(0).id()); 
	    Lima::Noeud n1 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(1).id()); 
	    Lima::Noeud n2 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(2).id()); 
	    Lima::Noeud n3 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(3).id()); 
	    Lima::Noeud n4 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(4).id()); 
	    Lima::Noeud n5 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(5).id()); 
	    Lima::Polyedre poly (maya.polyedre(ElementMinMax).id(),n0,n1,n2,n3,n4,n5);
	    mayaMinMax[iminmax].ajouter(poly);
	    mayaMinMax[iminmax].ajouter(vol);
	    vol.ajouter(poly);
	  }
	  break;
	case 8:
	  {
	    Lima::Noeud n0 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(0).id()); 
	    Lima::Noeud n1 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(1).id()); 
	    Lima::Noeud n2 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(2).id()); 
	    Lima::Noeud n3 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(3).id()); 
	    Lima::Noeud n4 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(4).id()); 
	    Lima::Noeud n5 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(5).id()); 
	    Lima::Noeud n6 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(6).id()); 
	    Lima::Noeud n7 = mayaMinMax[iminmax].noeud_id(maya.polyedre(ElementMinMax).noeud(7).id()); 
	    Lima::Polyedre poly (maya.polyedre(ElementMinMax).id(),n0,n1,n2,n3,n4,n5,n6,n7);
	    mayaMinMax[iminmax].ajouter(poly);
	    mayaMinMax[iminmax].ajouter(vol);
	    vol.ajouter(poly);
	  }
	  break;
	default:
	  break;
	}
      mayaMinMax[iminmax].ecrire(fichier); 
    }
}
