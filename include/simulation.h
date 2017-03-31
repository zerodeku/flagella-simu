//
// Created by zerodeku on 7/8/2016.
//

#ifndef FLAGELLA_SIMU_SIMULATION_H
#define FLAGELLA_SIMU_SIMULATION_H

namespace flagella {

struct Param {
  double l;
  double ds;
  double Rh;
  double beta;
  double gamma;
  double np;
  double cross_dist;
  double ka;
  double kb;
  double kc;
  double bodyK;
  double flagella_dist;
};

class Simulation {
  public:
    virtual void run() = 0;
    virtual void run_parallel(int n_threads) = 0;
};
}

#endif //FLAGELLA_SIMU_SIMULATION_H
