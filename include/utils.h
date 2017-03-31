//
// Created by zerodeku on 7/3/2016.
//

#ifndef FLAGELLA_SIMU_UTILS_H
#define FLAGELLA_SIMU_UTILS_H

#include <math.h>
#include <armadillo>

using namespace arma;

namespace flagella {

class Utils {
  public:
    static double const PI;

    static void quat_rotate(double const *vin, double const *axis,
                            double theta, double *vout);

    static double find_dist(double const *u, double const *v) {
      double tmp = 0;
      for (int i = 0; i < 3; i++) {
        tmp += pow(u[i] - v[i], 2);
      }
      return sqrt(tmp);
    }

    static void cross_prod(double const *u, double const *v,
                           double *vout) {
      vout[0] = u[1] * v[2] - u[2] * v[1];
      vout[1] = u[2] * v[0] - u[0] * v[2];
      vout[2] = u[0] * v[1] - u[1] * v[0];
    }
    static double dot_prod(double const *u, double const *v) {
      return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    }

    static void vect_diff(double const *u, double const *v,
                           double *vout) {
      vout[0] = u[0] - v[0];
      vout[1] = u[1] - v[1];
      vout[2] = u[2] - v[2];
    }

    static void add_vect(double *u, double const *v) {
      u[0] += v[0];
      u[1] += v[1];
      u[2] += v[2];
    }

    static double norm(double const *u) {
      return sqrt(pow(u[0], 2) + pow(u[1], 2) +pow(u[2], 2));
    }

    static void normalize(double *u) {
      double magn = norm(u);
      u[0] /= magn;
      u[1] /= magn;
      u[2] /= magn;
    }

    static void scale(double *u, double coeff) {
      u[0] *= coeff;
      u[1] *= coeff;
      u[2] *= coeff;
    }

    static void scale(double const *u, double coeff,
                      double *vout) {
      vout[0] = u[0] * coeff;
      vout[1] = u[1] * coeff;
      vout[2] = u[2] * coeff;
    }

    static void show_row(const rowvec &v) {
      cout << v(0) << ", " << v(1) << ", " << v(2) << endl;
    }

    static void show_col(const vec &v) {
      for (int i = 0; i < v.n_rows; i++) {
        cout << v(i) << endl;
      }
    }

    static void show_mat(const mat &m) {
      for (int i = 0; i < m.n_rows; i++) {
        show_row(m.row(i));
      }
    }

    static void copy_array(double *vin, double *vout) {
      for (int i = 0; i < 3; i++) {
        vout[i] = vin[i];
      }
    }

};


}

#endif //FLAGELLA_SIMU_UTILS_H
