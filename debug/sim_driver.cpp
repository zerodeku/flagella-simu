#include <iostream>
#include <ctime>
#include <windows.h>

#include <simulation.h>
#include <single_flagella.h>
#include <multi_flagella.h>

using namespace std;
using namespace flagella;

int main(int argc, char *argv[]) {
  int n_iters = 1;
  int n_threads = 4;
  bool parallel = false;
  bool single_simu = true;
  if (argc >= 2) {
    n_iters = atoi(argv[1]);
    if (argc == 3) {
      n_threads = atoi(argv[2]);
      parallel = true;
    }
  }
  cout << "Perform " << n_iters << " iterations." << endl;
  clock_t start = clock();
  Simulation *sim;
  if (single_simu) {
    sim = new SingleFlagella(n_iters);
  } else {
    sim = new MultiFlagella(n_iters);
  }
//  if (!parallel) {
//    sim->run();
//  } else {
//    sim->run_parallel(n_threads);
//  }
  sim->run_parallel(n_threads);

  double duration;
  duration = (clock() - start) / (double) CLOCKS_PER_SEC;
  cout << "Elapsed time: " << duration << "s" << endl;
  char buf[100];
  GetCurrentDirectory(100, buf);
  cout << "Result has been writen to folder: " << buf << endl;
  return 0;
}