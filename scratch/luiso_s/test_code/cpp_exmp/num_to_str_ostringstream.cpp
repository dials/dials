#include <sstream>     // for ostringstream
#include <string>
#include <iostream>

int main()
{
  const float fltValue = 1234.567890f;
  std::ostringstream ostr;

  ostr << fltValue;

  std::string strValue = ostr.str();

  std::cout << "fltValue = " << fltValue << std::endl;
  std::cout << "strValue = " << strValue << std::endl;

  return 0;
}
