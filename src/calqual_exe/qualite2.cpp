// Application de test de qualite de maillage

#include <fstream>
#include <cmath>

extern "C"
{
#include <unistd.h>
}

#include "Lima/erreur.h"
#include "Lima/polyedre.h"
#include "Lima/volume.h"
#include "Qualif.h"

//#include "mathimf.h"
std::string
  type_polyedre[] =
{
"tetraedre", "pyramide", "prisme", "hexaedre", "heptaedre", "octoedre"};


void
analyse (Lima::Volume volume, int numCritere)
{
  std::cout << "Analyse de qualite pour le volume " << volume.nom ()
    << std::endl;

  int
    nbElts =
    volume.
    nb_polyedres ();
  int
    nbEltsMesures =
    nbElts;
  double *
    Resultats =
    new double[nbElts];

  Qualif::AnalyseVolume (volume, Qualif::Critere (numCritere), Resultats);


  double sx[Lima::Polyedre::maxPolyedreType];
  double sx2[Lima::Polyedre::maxPolyedreType];
  int nbr[Lima::Polyedre::maxPolyedreType];
  int nbval[Lima::Polyedre::maxPolyedreType];
  int val_min[Lima::Polyedre::maxPolyedreType];
  int val_max[Lima::Polyedre::maxPolyedreType];

  for (int i = 0; i < Lima::Polyedre::maxPolyedreType; ++i)
    {
      sx[i] = sx2[i] = 0.0;
      nbr[i] = nbval[i] = 0;
    }

  double total_sx, total_sx2;
  int total_nbr, total_nbval;

  total_sx = total_sx2 = 0.0;
  total_nbr = total_nbval = 0;

  for (int pol_id = 0; pol_id < nbElts; pol_id++)
    {
      double
	res =
	Resultats[pol_id];
      Lima::Polyedre::PolyedreType poly_type =
	volume.polyedre (pol_id).type ();
      bool valid = (std::isnan (res) == false);

      nbr[poly_type]++;
      total_nbr++;

      if (valid)
	{
	  sx[poly_type] += res;
	  sx2[poly_type] += res * res;

	  total_sx += res;
	  total_sx2 += res * res;
	  nbval[poly_type]++;
	  if (nbval[poly_type] == 1)
	    {
	      val_min[poly_type] = val_max[poly_type] = pol_id;
	    }
	  else if (res < Resultats[val_min[poly_type]])
	    val_min[poly_type] = pol_id;
	  else if (res > Resultats[val_max[poly_type]])
	    val_max[poly_type] = pol_id;
	}
    }

  for (int type = 0; type < Lima::Polyedre::maxPolyedreType; type++)
    {
      if (nbr[type] != 0)
	{
	  std::cout << "Resultats pour " << type_polyedre[type] << std::endl;
	  std::cout << "nbr = " << nbr[type] << "  valides = " << nbval[type]
	    << std::endl;
	  if (nbval[type] != 0)
	    {
	      std::cout << "moyenne " << sx[type] / (double) nbval[type];
	      std::cout << " ecart type ";
	      std::cout << std::
		sqrt ((sx2[type] -
		       sx[type] * sx[type] / nbval[type]) / nbval[type]);


	      std::cout << std::endl;

	      std::cout << "min = " << Resultats[val_min[type]];
	      std::cout << " ( " << volume.polyedre (val_min[type]).id () << ")";
	      std::cout << " max = " << Resultats[val_max[type]];
	      std::cout << " ( " << volume.polyedre (val_max[type]).id () << ")";
	      std::cout << std::endl;
	      std::cout << std::endl;

	    }
	}
    }
  delete[]Resultats;
}

int
main (int argc, char *argv[])
{

  int
    numCritere =  0;
  int
    opt =::getopt (argc, argv, "c:");
  if (opt == 'c')
    numCritere =::atoi (optarg);

  std::cout << "     ************************************" << std::endl;
  std::cout << "         MESURE DE QUALITE DE MAILLAGE   " << std::endl;
  std::cout << "         Bibliotheque Qualif " << Qualif::QualifVersion () << std::endl;
  std::cout << "     ************************************" << std::endl;
  std::cout << " " << std::endl;


  if (optind == argc)
    {
      std::cerr << "Usage "
	<< argv[0] << " [-c critere] nom_maillage (volume)* " << std::endl;
      exit (-1);
    }

  //======================================
  //         LECTURE DU MAILLAGE
  //======================================
  Lima::Maillage maya;
  char  rep;

  try
  {
    maya.lire (argv[optind++]);

    if (numCritere == 0)
      {
	Qualif::AfficheCriteres ();
	std::cout << " -> NUMERO DU CRITERE A APPLIQUER ? " << std::endl;

	std::cin >> numCritere;
	std::cout << " " << std::endl;
      }
    std::cout << " APPLICATION DU CRITERE : " << Qualif::
      CRITERESTR[numCritere] << std::endl;

    if (optind < argc)
      {
	for (int i = optind; i < argc; i++)
	  {
	    Lima::Volume volume = maya.volume (argv[i]);
	    analyse (volume, numCritere);
	  }
      }
    else
      {
	for (int vol = 0; vol < maya.nb_volumes (); vol++)
	  {
	    Lima::Volume volume = maya.volume (vol);
	    analyse (volume, numCritere);
	  }
      }
  }
  catch (Lima::erreur & exc)
  {
    std::cerr << "Lima reporte l'erreur suivante : "
      << exc.what () << std::endl;
    exit (-1);
  }
}
