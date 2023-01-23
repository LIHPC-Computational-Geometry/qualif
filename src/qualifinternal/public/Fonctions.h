#ifndef DEF_FONCTIONS_H
#define DEF_FONCTIONS_H

namespace Qualif {

//========================================================
//              FONCTION RENVOYANT L'ADRESSE 
//           DU PLUS GRAND ELEMENT D'UN TABLEAU
//                  OU LE POINTEUR NUL
//========================================================
inline double* PlusGrandElement(double *tab, int n) 
{
  double *ptr_max = 0;
  double Element_max = -DBL_MAX;
  for (int i =0; i < n; i++)
    if (tab[i] >= Element_max) 
      {
	ptr_max     = &tab[i];
	Element_max = tab[i];
      } 
  return ptr_max;
}

//========================================================
//              FONCTION RENVOYANT L'ADRESSE 
//           DU PLUS PETIT ELEMENT D'UN TABLEAU
//                  OU LE POINTEUR NUL
//========================================================
inline double* PlusPetitElement(double *tab, int n) 
{
  double *ptr_min = 0;
  double Element_min = DBL_MAX;
  for (int i =0; i < n; i++)
    if (tab[i] <= Element_min) 
      {
	ptr_min     = &tab[i];
	Element_min = tab[i];
      } 
  return ptr_min;
}


//========================================================
//          MAX DE DEUX NOMBRES
//          MIN DE DEUX NOMBRES
//========================================================

inline double min (double x, double y)
{       
       if (x < y)
        return x;
       return y;
}

inline double max (double x, double y)
{       
       if (x < y)
        return y;
       return x;
}

}

#endif

