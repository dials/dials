/* atof example: sine calculator */
#include <stdio.h>      /* printf, fgets */
#include <stdlib.h>     /* atof */
#include <math.h>       /* sin */

int main ()
{

  double n,m;
  double pi=3.1415926535;
  char buffer[256];

  printf ("Enter degrees: ");
  fgets (buffer,256,stdin);
  n = atof (buffer);
  m = sin (n*pi/180);
  printf ("The sine of %f degrees is %f\n" , n, m);

  return 0;
}
