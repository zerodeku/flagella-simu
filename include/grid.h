//
// Created by zerodeku on 7/2/2016.
//

#ifndef FLAGELLA_SIMU_GRID_H
#define FLAGELLA_SIMU_GRID_H

#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

namespace flagella {

class Grid {
  public:
    Grid() { }
    Grid(double x, double y, double z) {
      pos_ = new double[3]{x, y, z};
      adj_idx_ = new vector<int>();
      stiff_const_ = new vector<double>();
      rest_len_ = new vector<double>();
      adj_idx_->clear();
      stiff_const_->clear();
      rest_len_->clear();
//      std::cout << "Grid created!" << std::endl;
    }

    ~Grid() {
      delete pos_;
      delete adj_idx_;
      delete stiff_const_;
      delete rest_len_;
//      std::cout << "Grid destroyed!" << std::endl;
    }

    double *pos_;
    vector<int> *adj_idx_;
    vector<double> *stiff_const_;
    vector<double> *rest_len_;

    // Initialize the grid conditions
    static void BuildGrids(struct Param p, Grid **grids, int n);
    // Update implementation, making the configuration centrosymmetric
    static void BuildGrids1(struct Param p, Grid **grids, int n);
    // Build grids for multiple flagella
    static void BuildGridsMultipleFlagella(struct Param p, Grid **grids, int n,
                                           int n_flagella, bool rotate_flag);
};

}

#endif //FLAGELLA_SIMU_GRID_H
