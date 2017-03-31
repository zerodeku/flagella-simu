//
// Created by zerodeku on 7/3/2016.
//

#include <utils.h>

using namespace std;

namespace flagella {

double const Utils::PI = 3.1415926535897;

void Utils::quat_rotate(double const *vin, double const *axis,
                        double theta, double *vout) {
  int len = 4;
  double q[len];
  q[0] = (double)cos(theta / 2);
  double norm = (double)sqrt(pow(axis[0], 2) + pow(axis[1], 2) +pow(axis[2], 2));
  for (int i = 1; i < len; i++) {
    q[i] = (double)(axis[i - 1] * sin(theta / 2) / norm);
  }

  double rot_mat[3][3] = {
      {(double)(1-2*pow(q[2],2)-2*pow(q[3],2)), 2*(q[1]*q[2]-q[0]*q[3]), 2*(q[1]*q[3]+q[0]*q[2])},
      {2*(q[1]*q[2]+q[0]*q[3]), (double)(1-2*pow(q[1],2)-2*pow(q[3],2)), 2*(q[2]*q[3]-q[0]*q[1])},
      {2*(q[1]*q[3]-q[0]*q[2]), 2*(q[2]*q[3]+q[0]*q[1]), (double)(1-2*pow(q[1],2)-2*pow(q[2],2))}
  };

  for (int i = 0; i < 3; i++) {
    vout[i] = 0;
    for (int j = 0; j < 3; j++) {
      vout[i] += rot_mat[i][j] * vin[j];
    }
  }

  // normalize the output vector
  norm = (double)sqrt(pow(vout[0], 2) + pow(vout[1], 2) +pow(vout[2], 2));
  for (int i = 0; i < 3; i++) {
    vout[i] /= norm;
  }
}

}