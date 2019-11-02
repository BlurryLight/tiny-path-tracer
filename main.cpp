#include "vec3.h"
#include <iostream>

using namespace std;

int main()
{
  int nx = 200;
  int ny = 100;
  std::cout << "P3\n" << nx << " " << ny << "\n255\n";
  for (int j = ny - 1; j >= 0; j--) {
    for (int i = 0; i < nx; i++) {
      // left-top (0.0,0.9,0.2)//green
      // left-b (0.0,0.0,0.2)//blue
      // right-t (1.0,1,0.2)//yellow

      vec3 color(float(i) / float(nx), float(j) / float(ny), 0.2f);

      int red = int(255.99f * color.r());
      int green = int(255.99f * color.g());
      int blue = int(255.99f * color.b());
      std::cout << red << ' ' << green << ' ' << blue << ' ' << '\n';
    }
  }
  return 0;
}
