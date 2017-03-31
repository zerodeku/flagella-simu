//
// Created by zerodeku on 7/3/2016.
//

#ifndef FLAGELLA_SIMU_SINGLE_FLAGELLA_H
#define FLAGELLA_SIMU_SINGLE_FLAGELLA_H

#include <simulation.h>
#include <grid.h>
#include <armadillo>
#include <math.h>
#include <fstream>
#include <boost/thread/thread.hpp>

using namespace arma;

namespace flagella {

class SingleFlagella: public Simulation {
  public:
    int num_iters_ = 0;
    int iter_ = 0;
    Param params_ = {
        1.3,
        0.025,
        0.04,
        10,
        0.1,
        3,
        0.026,
        12,
        12,
        12,
        12,
        0.09,
    };
    Grid **grids_;
    int n_ = round(params_.l / params_.ds) + 1;
    int num_grids_ = n_ * 3;
    double dt_ = 0.001;
    double deltaS_ = 0.039;
    double deltaR_ = 0.052;
    double motorT_ = 0.002;
    mat pos_;
    mat f_;
    mat pos_next_;
    rowvec posL_ = rowvec(3);
    rowvec vectL_ = rowvec(3);
    vector<mat> all_pos_;
    fstream pos_out_;
    fstream f_out_;
    int n_threads_ = 4;

    // constructor
    SingleFlagella(int num_iters) {
      num_iters_ = num_iters;
      pos_ = mat(num_grids_, 3);
      pos_next_ = mat(num_grids_, 3);
      f_ = mat(num_grids_, 3);
      pos_.zeros();
      f_.zeros();
      pos_out_.open("pos_result.txt", fstream::out);
      f_out_.open("f_result.txt", fstream::out);
    }
    void build_grids();
    void find_next_force();
    void find_next_pos();
    void find_next_pos_parallel();
    void find_next_force_parallel();
    void update_torque();
    void update_torque1();
    void emit_result();
    void find_pos_in_range(int start, int end);
    void find_force_in_range(int start, int end);

    // run the simulation
    void run();
    void run_parallel(int n_threads);
};

}

#endif //FLAGELLA_SIMU_SINGE_FLAGELLA_H
