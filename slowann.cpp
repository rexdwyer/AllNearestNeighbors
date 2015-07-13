#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#define REALBIG 9e99  
using namespace std;


class Site { public: int index; double x, y;};

main() {

  vector<Site*> sites;
  std::string line;
  int ii = 1;
  while (getline(cin,line)) {
    istringstream iss(line);
    Site * site = new Site;
    site->index = ii++;
    iss >> site->x >> site->y;
    sites.push_back(site);
  }

  for (vector<Site*>::iterator i = sites.begin(); i != sites.end(); ++i) {
    Site *sitei = *i;
    int bestj;
    double best = 9e99;
    for (vector<Site*>::iterator j = sites.begin(); j != sites.end(); ++j) {
      Site *sitej = *j;
      if (sitei->index == sitej->index) continue;
      double xdiff = sitei->x - sitej->x;
      double ydiff = sitei->y - sitej->y;
      double dist2 = xdiff*xdiff + ydiff*ydiff;
      if (dist2 < best) { best = dist2; bestj = sitej->index; }
    }
    cout << sitei->index << " " << bestj << endl;
  }
}
