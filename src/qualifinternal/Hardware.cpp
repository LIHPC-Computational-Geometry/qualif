//
// Fichier pour generer les bons resultats sur NAN et Infini
//

#include "IQualif.h"

namespace Qualif {
 
#ifndef WORD_LENGTH_64
static const unsigned long signbit_mask          = 0x80000000;
static const unsigned long clear_signbit_mask    = 0x7fffffff;
static const unsigned long exponent_mask         = 0x7ff00000;
static const unsigned long single_precision_mask = 0x07ffffff;
static const unsigned long high_mantissa_mask    = 0x000fffff;
#else
static const unsigned long signbit_mask          = 0x8000000000000000;
static const unsigned long clear_signbit_mask    = 0x7fffffffffffffff;
static const unsigned long exponent_mask         = 0x7ff0000000000000;
static const unsigned long single_precision_mask = 0x0000000007ffffff;
static const unsigned long mantissa_mask         = 0x000fffffffffffff;
#endif

#ifndef WORD_LENGTH_64

inline double compose_parts(int sign_1, unsigned long exp_11,
                     unsigned long most_sig_20, unsigned long least_sig_32)
{
  double a;

  unsigned long high_32=0;
  if (sign_1) 
    high_32 |= signbit_mask;
  high_32 |= (exp_11 << 20);
  high_32 |= most_sig_20;
  
  unsigned long *p;
  p=(unsigned long*)&a;
#ifndef LITTLE_ENDIAN_MACHINE
  *p =high_32; 
   p++; 
  *p =least_sig_32;
#else
  *p =least_sig_32; 
   p++; 
  *p =high_32;
#endif

  return a;
}


void read_parts(const double& a,
                int& sign_1, long& exp_11,
                unsigned long& most_sig_20, unsigned long& least_sig_32)
{
  unsigned long *p = (unsigned long*) &a;
  unsigned long high_32;
#ifndef LITTLE_ENDIAN_MACHINE
  high_32 = *p; 
  p++; 
  least_sig_32 = *p;
#else
  least_sig_32 = *p; 
  p++; 
  high_32 = *p;
#endif

  sign_1 = (int) (high_32 & signbit_mask) >> 31;
  exp_11 = (high_32 & exponent_mask) >> 20;
  most_sig_20 = (high_32 & 0x000fffff);
}

long binary_equal(const double& x, const double& y)
{
  unsigned long* x_ptr = (unsigned long*) &x;
  unsigned long* y_ptr = (unsigned long*) &y;
  if (*x_ptr != *y_ptr)
    return 0;
  x_ptr++; y_ptr++;
  if (*x_ptr != *y_ptr)
    return 0;
  return 1;
}

#else

double compose_parts(int sign_1, unsigned long exp_11,
                     unsigned long most_sig_20, unsigned long least_sig_32)
{
  double a;

  unsigned long word_64 = (exp_11 << 20) + most_sig_20;
  word_64 = (word_64 << 32) + least_sig_32;
  if (sign_1) 
    word_64 |= signbit_mask;
  
  unsigned long *p;
  p=(unsigned long*)&a;
  *p = word_64;

  return a;
}


void read_parts(const double& a,
                int& sign_1, long& exp_11,
                unsigned long& most_sig_20, unsigned long& least_sig_32)
{
  unsigned long *p = (unsigned long*) &a;
  unsigned long word_64 = *p;
  sign_1 = (int) (word_64 >> 63);
  exp_11 = (word_64 & exponent_mask) >> 52;
  most_sig_20  = (word_64 & 0x000fffff00000000) >> 32;
  least_sig_32 = word_64 & 0x00000000ffffffff;
}

long binary_equal(const double& x, const double& y)
{
  unsigned long* x_ptr = (unsigned long*) &x;
  unsigned long* y_ptr = (unsigned long*) &y;
  return (*x_ptr == *y_ptr);
}

#endif



const double double_min=compose_parts(0,1,0,0);

const double Not_A_Number=compose_parts(0,2047,0,1);
const double Infini=compose_parts(0,2047,0,0);

const double pZero_double=compose_parts(0,0,0,0);
const double nZero_double=compose_parts(1,0,0,0);



enum hardware_type { B_ENDIAN=0, L_ENDIAN=1, LENGTH_64 = 2};

hardware_type check_status()
{
   int size = sizeof(double);

#ifndef WORD_LENGTH_64
   if (size == sizeof(long))
   {
      std::cout << 
        "error in compilation of hardware.cc: use flag -DWORD_LENGTH_64\n";
      exit(-1);
   }
   const static double x = 1;
   long *px = (long *) &x;
   
#ifndef LITTLE_ENDIAN_MACHINE
   if (*(++px))
   {
     std::cout <<   
       "error in compilation of hardware.cc: use flag -DLITTLE_ENDIAN_MACHINE\n";
     exit(-1);
   }
   return B_ENDIAN;
#else
   if (*px)
   {
     std::cout << 
        "error in compilation of hardware.cc: don't use flag -DLITTLE_ENDIAN_MACHINE\n";
     exit(-1);
   }
   return L_ENDIAN;
#endif


#else
   if (sizeof(long) < size)
   {
      std::cout << 
        "error in compilation of hardware.cc: don't use flag -DWORD_LENGTH_64\n";
      exit(-1);
   }
   return LENGTH_64;
#endif

}

hardware_type type = check_status();

}
