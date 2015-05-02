/* sprintf example */
#include <stdio.h>
#include <iostream>
//#include <string>     // std::string, std::to_string
int main ()
{
  char asc_str [50];

  int i_nm = 1235434567;
  sprintf (asc_str, "%d ", i_nm);
  std::cout << "\nasc_str = " << asc_str << "\n\n";

  double dl_nm = 8876576.465654654;
  sprintf (asc_str, "%f ", dl_nm);
  std::cout << "\nasc_str = " << asc_str << "\n\n";

  return 0;
}
