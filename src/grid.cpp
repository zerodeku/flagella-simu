//
// Created by zerodeku on 7/2/2016.
//

#include <iostream>

#include <grid.h>
#include <utils.h>
#include <single_flagella.h>

namespace flagella {

void Grid::BuildGrids(struct Param p, Grid **grids, int n) {
  const double pi = Utils::PI;
  double s[n], Rs, Rsp, Xs, Xsp, Ys, Ysp, Zs, Zsp, pre_Xs;
  double Xs1, Ys1, Zs1, Xs2, Ys2, Zs2, mgn, nXs1, nZs1;

  s[0] = 0;
  for (int i = 1; i < n; i++) {
    s[i] = s[i - 1] + p.ds;
  }
  Xs = 0;
  pre_Xs = 0;
  Xsp = 0;
  for (int i = 0; i < n; i++) {
    Rs = p.Rh * (1 / pi * atan(p.beta * (s[i] / p.l - p.gamma)) + 0.5);
    Rsp = p.Rh / pi / (1 + pow((p.beta * (s[i] / p.l - p.gamma)), 2)) *
          p.beta / p.l;
    Ys = -Rs * cos(2 * pi * p.np * s[i] / p.l);
    Ysp = -Rsp * cos(2 * pi * p.np * s[i] / p.l) + Rs *
               sin(2 * pi * p.np * s[i] / p.l) * 2 * pi * p.np / p.l;
    Zs = Rs * sin(2 * pi * p.np * s[i] / p.l);
    Zsp = Rsp * sin(2 * pi * p.np * s[i] / p.l) + Rs *
              cos(2 * pi * p.np * s[i] / p.l) * 2 * pi * p.np / p.l;

    if (i > 0) {
      Xs = pre_Xs + Xsp * p.ds;
      pre_Xs = Xs;
    }
    Xsp = sqrt(1 - pow(Ysp, 2) - pow(Zsp, 2));


    mgn = sqrt(pow(Xsp, 2) + pow(Zsp, 2));
    nXs1 = Zsp / mgn;
    nZs1 = -Xsp / mgn;
    Xs1 = Xs + nXs1 * p.cross_dist;
    Ys1 = Ys;
    Zs1 = Zs + nZs1 * p.cross_dist;

    double vin[3] = {nXs1, 0, nZs1};
    double vout[3];
    double axis[3] = {Xsp, Ysp, Zsp};
    Utils::quat_rotate(vin, axis, pi / 3, vout);
    Xs2 = vout[0] * p.cross_dist + Xs;
    Ys2 = vout[1] * p.cross_dist + Ys;
    Zs2 = vout[2] * p.cross_dist + Zs;

    grids[i * 3] = new Grid(Xs, Ys, Zs);
    grids[i * 3 + 1] = new Grid(Xs1, Ys1, Zs1);
    grids[i * 3 + 2] = new Grid(Xs2, Ys2, Zs2);
  }

  int len = n * 3;
  for (int i = 0; i < len; i += 3) {
    grids[i]->adj_idx_->push_back(i + 1);
    grids[i]->adj_idx_->push_back(i + 2);
    grids[i]->stiff_const_->push_back(p.ka);
    grids[i]->stiff_const_->push_back(p.ka);
    if (i - 3 >= 0) {
      grids[i]->adj_idx_->push_back(i - 1);
      grids[i]->adj_idx_->push_back(i - 2);
      grids[i]->adj_idx_->push_back(i - 3);

      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kb);
    }
    if (i + 5 < len) {
      grids[i]->adj_idx_->push_back(i + 3);
      grids[i]->adj_idx_->push_back(i + 4);
      grids[i]->adj_idx_->push_back(i + 5);

      grids[i]->stiff_const_->push_back(p.kb);
      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kc);
    }
  }
  for (int i = 1; i < len; i += 3) {
    grids[i]->adj_idx_->push_back(i - 1);
    grids[i]->adj_idx_->push_back(i + 1);
    grids[i]->stiff_const_->push_back(p.ka);
    grids[i]->stiff_const_->push_back(p.ka);
    if (i - 4 >= 0) {
      grids[i]->adj_idx_->push_back(i - 2);
      grids[i]->adj_idx_->push_back(i - 3);
      grids[i]->adj_idx_->push_back(i - 4);

      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kb);
      grids[i]->stiff_const_->push_back(p.kc);
    }
    if (i + 4 < len) {
      grids[i]->adj_idx_->push_back(i + 2);
      grids[i]->adj_idx_->push_back(i + 3);
      grids[i]->adj_idx_->push_back(i + 4);

      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kb);
      grids[i]->stiff_const_->push_back(p.kc);
    }
  }
  for (int i = 2; i < len; i += 3) {
    grids[i]->adj_idx_->push_back(i - 1);
    grids[i]->adj_idx_->push_back(i - 2);
    grids[i]->stiff_const_->push_back(p.ka);
    grids[i]->stiff_const_->push_back(p.ka);
    if (i - 5 >= 0) {
      grids[i]->adj_idx_->push_back(i - 3);
      grids[i]->adj_idx_->push_back(i - 4);
      grids[i]->adj_idx_->push_back(i - 5);

      grids[i]->stiff_const_->push_back(p.kb);
      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kc);
    }
    if (i + 3 < len) {
      grids[i]->adj_idx_->push_back(i + 1);
      grids[i]->adj_idx_->push_back(i + 2);
      grids[i]->adj_idx_->push_back(i + 3);

      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kb);
    }
  }

  for (int i = 0; i < len; i++) {
    for (int idx : *(grids[i]->adj_idx_)) {
      grids[i]->rest_len_->push_back(Utils::find_dist(grids[i]->pos_,
                                                     grids[idx]->pos_));
    }
  }
}

void Grid::BuildGrids1(struct Param p, Grid **grids, int n) {
  const double pi = Utils::PI;
  double s[n], Rs, Rsp, Xs, Xsp, Ys, Ysp, Zs, Zsp, pre_Xs;
  double Xs1, Ys1, Zs1, Xs2, Ys2, Zs2, mgn, nXs1, nZs1;

  s[0] = 0;
  for (int i = 1; i < n; i++) {
    s[i] = s[i - 1] + p.ds;
  }
  Xs = 0;
  pre_Xs = 0;
  Xsp = 0;
  for (int i = 0; i < n; i++) {
    Rs = p.Rh * (1 / pi * atan(p.beta * (s[i] / p.l - p.gamma)) + 0.5);
    Rsp = p.Rh / pi / (1 + pow((p.beta * (s[i] / p.l - p.gamma)), 2)) *
          p.beta / p.l;
    Ys = -Rs * cos(2 * pi * p.np * s[i] / p.l);
    Ysp = -Rsp * cos(2 * pi * p.np * s[i] / p.l) + Rs *
                                                   sin(2 * pi * p.np * s[i] / p.l) * 2 * pi * p.np / p.l;
    Zs = Rs * sin(2 * pi * p.np * s[i] / p.l);
    Zsp = Rsp * sin(2 * pi * p.np * s[i] / p.l) + Rs *
                                                  cos(2 * pi * p.np * s[i] / p.l) * 2 * pi * p.np / p.l;

    if (i > 0) {
      Xs = pre_Xs + Xsp * p.ds;
      pre_Xs = Xs;
    }
    Xsp = sqrt(1 - pow(Ysp, 2) - pow(Zsp, 2));


    double x = Ys * Ysp + Zs * Zsp;
    double y = -Ys;
    double z = -Zs;
    double vin[3] = {-x, -y, -z};
    Utils::normalize(vin);
    double vout[3];
    double axis[3] = {Xsp, Ysp, Zsp};

    Utils::quat_rotate(vin, axis, - pi / 6, vout);
//    Utils::normalize(vout);
    Xs1 = vout[0] * p.cross_dist + Xs;
    Ys1 = vout[1] * p.cross_dist + Ys;
    Zs1 = vout[2] * p.cross_dist + Zs;

    Utils::quat_rotate(vin, axis, pi / 6, vout);
//    Utils::normalize(vout);
    Xs2 = vout[0] * p.cross_dist + Xs;
    Ys2 = vout[1] * p.cross_dist + Ys;
    Zs2 = vout[2] * p.cross_dist + Zs;

    grids[i * 3] = new Grid(Xs, Ys, Zs);
    grids[i * 3 + 1] = new Grid(Xs1, Ys1, Zs1);
    grids[i * 3 + 2] = new Grid(Xs2, Ys2, Zs2);
  }

  int len = n * 3;
  for (int i = 0; i < len; i += 3) {
    grids[i]->adj_idx_->push_back(i + 1);
    grids[i]->adj_idx_->push_back(i + 2);
    grids[i]->stiff_const_->push_back(p.ka);
    grids[i]->stiff_const_->push_back(p.ka);
    if (i - 3 >= 0) {
      grids[i]->adj_idx_->push_back(i - 1);
      grids[i]->adj_idx_->push_back(i - 2);
      grids[i]->adj_idx_->push_back(i - 3);

      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kb);
    }
    if (i + 5 < len) {
      grids[i]->adj_idx_->push_back(i + 3);
      grids[i]->adj_idx_->push_back(i + 4);
      grids[i]->adj_idx_->push_back(i + 5);

      grids[i]->stiff_const_->push_back(p.kb);
      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kc);
    }
  }
  for (int i = 1; i < len; i += 3) {
    grids[i]->adj_idx_->push_back(i - 1);
    grids[i]->adj_idx_->push_back(i + 1);
    grids[i]->stiff_const_->push_back(p.ka);
    grids[i]->stiff_const_->push_back(p.ka);
    if (i - 4 >= 0) {
      grids[i]->adj_idx_->push_back(i - 2);
      grids[i]->adj_idx_->push_back(i - 3);
      grids[i]->adj_idx_->push_back(i - 4);

      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kb);
      grids[i]->stiff_const_->push_back(p.kc);
    }
    if (i + 4 < len) {
      grids[i]->adj_idx_->push_back(i + 2);
      grids[i]->adj_idx_->push_back(i + 3);
      grids[i]->adj_idx_->push_back(i + 4);

      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kb);
      grids[i]->stiff_const_->push_back(p.kc);
    }
  }
  for (int i = 2; i < len; i += 3) {
    grids[i]->adj_idx_->push_back(i - 1);
    grids[i]->adj_idx_->push_back(i - 2);
    grids[i]->stiff_const_->push_back(p.ka);
    grids[i]->stiff_const_->push_back(p.ka);
    if (i - 5 >= 0) {
      grids[i]->adj_idx_->push_back(i - 3);
      grids[i]->adj_idx_->push_back(i - 4);
      grids[i]->adj_idx_->push_back(i - 5);

      grids[i]->stiff_const_->push_back(p.kb);
      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kc);
    }
    if (i + 3 < len) {
      grids[i]->adj_idx_->push_back(i + 1);
      grids[i]->adj_idx_->push_back(i + 2);
      grids[i]->adj_idx_->push_back(i + 3);

      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kc);
      grids[i]->stiff_const_->push_back(p.kb);
    }
  }

  for (int i = 0; i < len; i++) {
    for (int idx : *(grids[i]->adj_idx_)) {
      grids[i]->rest_len_->push_back(Utils::find_dist(grids[i]->pos_,
                                                      grids[idx]->pos_));
    }
  }
}

void Grid::BuildGridsMultipleFlagella(struct Param p, Grid **grids, int n,
                                      int n_flagella, bool rotate_flag) {
  if (n_flagella == 2) {
    BuildGrids1(p, grids, n);
    BuildGrids1(p, &grids[3 * n], n);
    for (int i = n * 3; i < n * 3 * 2; i++) {
      for (int j = 0; j < grids[i]->adj_idx_->size(); j++) {
        grids[i]->adj_idx_->at(j) = grids[i]->adj_idx_->at(j) + n * 3;
      }
      grids[i]->pos_[2] -= p.flagella_dist;
    }

    int idx_shift = 3 * n;
    double tmp[3] = {0, 0, 0};
    double base1[3];
    int base1_idx = n * 3 * 2, base2_idx = n * 3 * 2 + 1;
    Utils::add_vect(tmp, grids[0]->pos_);
    Utils::add_vect(tmp, grids[1]->pos_);
    Utils::add_vect(tmp, grids[2]->pos_);
    Utils::scale(tmp, 1.0 / 3, base1);
    grids[base1_idx] = new Grid(base1[0], base1[1], base1[2]);
    tmp[0] = 0;
    tmp[1] = 0;
    tmp[2] = 0;
    Utils::add_vect(tmp, grids[n * 3]->pos_);
    Utils::add_vect(tmp, grids[n * 3 + 1]->pos_);
    Utils::add_vect(tmp, grids[n * 3 + 2]->pos_);
    Utils::scale(tmp, 1.0 / 3, base1);
    grids[base2_idx] = new Grid(base1[0], base1[1], base1[2]);

    for (int i = 0; i < 3; i++) {
      grids[i]->adj_idx_->push_back(base1_idx);
      grids[i]->stiff_const_->push_back(p.ka);
      grids[i]->rest_len_->push_back(Utils::find_dist(grids[i]->pos_,
                                                      grids[base1_idx]->pos_));
      grids[base1_idx]->adj_idx_->push_back(i);
      grids[base1_idx]->stiff_const_->push_back(p.ka);
      grids[base1_idx]->rest_len_->push_back(Utils::find_dist(grids[i]->pos_,
                                                      grids[base1_idx]->pos_));
      int j = i + n * 3;
      grids[j]->adj_idx_->push_back(base2_idx);
      grids[j]->stiff_const_->push_back(p.ka);
      grids[j]->rest_len_->push_back(Utils::find_dist(grids[j]->pos_,
                                                      grids[base2_idx]->pos_));
      grids[base2_idx]->adj_idx_->push_back(j);
      grids[base2_idx]->stiff_const_->push_back(p.ka);
      grids[base2_idx]->rest_len_->push_back(Utils::find_dist(grids[j]->pos_,
                                                              grids[base2_idx]->pos_));
    }

    grids[base1_idx]->adj_idx_->push_back(base2_idx);
    grids[base1_idx]->stiff_const_->push_back(p.bodyK);
    grids[base1_idx]->rest_len_->push_back(p.flagella_dist);
    grids[base2_idx]->adj_idx_->push_back(base1_idx);
    grids[base2_idx]->stiff_const_->push_back(p.bodyK);
    grids[base2_idx]->rest_len_->push_back(p.flagella_dist);

    if (rotate_flag) {
      double axis[3] = {0, 1, 0};
      double theta = Utils::PI / 3;
      double ref_point[3];
      double tmp_vect[3];
      double trans_vect[3];
      int ref_idx = n * 3 * 2 + 1;
      double length;
      Utils::copy_array(grids[ref_idx]->pos_, ref_point);
      length = Utils::norm(grids[ref_idx]->pos_);
      Utils::quat_rotate(ref_point, axis, theta, tmp_vect);
      Utils::scale(tmp_vect, length);
      Utils::vect_diff(ref_point, tmp_vect, trans_vect);

      int offset = n * 3;
      for (int i = offset; i < offset * 2; i++) {
        length = Utils::norm(grids[i]->pos_);
        Utils::quat_rotate(grids[i]->pos_, axis, theta, tmp_vect);
        Utils::scale(tmp_vect, length);
        Utils::add_vect(tmp_vect, trans_vect);
        Utils::copy_array(tmp_vect, grids[i]->pos_);
      }
    }
  }
}

}


