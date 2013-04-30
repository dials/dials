#include <iostream>
int main(){
  int data2d [10][10];
  for (int x=0;x<10;x++){
    for (int y=0;y<10;y++){
      std::cout << "x = "<<x<<"\ny = "<<y<<"\n";
      data2d[x][y]=x+y;
    }
  }
  std::cout<<data2d;

}
