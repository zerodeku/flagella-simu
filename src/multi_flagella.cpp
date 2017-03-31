//
// Created by zerodeku on 7/3/2016.
//

#include <multi_flagella.h>
#include <utils.h>

using namespace arma;

namespace flagella {

// Buid up grids of this iteration
void MultiFlagella::build_grids_double() {
  grids_ = new Grid*[num_grids_];
  Grid *g = new Grid(-1, -1, -1);
  for (int i = 0; i < num_grids_; i++) {
    grids_[i] = g;
  }
  Grid::BuildGridsMultipleFlagella(params_, grids_, n_, 2, rotate_flag_);
}

// find positionis of next iteration in parallel
void MultiFlagella::find_next_pos_parallel() {
  boost::thread_group find_pos_threads;
  int step = num_grids_ / n_threads_;
  int end;
  for (int i = 0; i < num_grids_; i += step) {
    end = (i + step > num_grids_) ? num_grids_ : i + step;
    find_pos_threads.create_thread(boost::bind(find_pos_in_range, this, i, end));
  }
  find_pos_threads.join_all();
  pos_ = pos_next_;
}

// find the next position in a range of grids
void MultiFlagella::find_pos_in_range(int start, int end) {
  for (int i = start; i < end; i++) {
    // rotlet induced vel
    rowvec v;
    rowvec vel = {0, 0, 0};
    for (int j = 0; j < posL_.n_rows; j++) {
      v = pos_.row(i) - posL_.row(j);
      double rr = pow(norm(v), 2);
      double coeff = (2 * rr + 5 * pow(deltaR_,2)) /
                     (16 * Utils::PI * pow(rr + pow(deltaR_,2), 2.5));
      vel += cross(vectL_.row(j), v) * coeff;
    }

    // stokeslet induced pos changes
    vec v0 = pos_.col(0) - pos_(i, 0);
    vec v1 = pos_.col(1) - pos_(i, 1);
    vec v2 = pos_.col(2) - pos_(i, 2);
    vec rs = pow(v0, 2) + pow(v1, 2) + pow(v2, 2);
    vec tmp = 8 * Utils::PI * pow(rs + pow(deltaS_,2), 1.5);
    vec coeffs = (rs + 2 * pow(deltaS_, 2)) / tmp;
    vel(0) += sum(f_.col(0) % coeffs);
    vel(1) += sum(f_.col(1) % coeffs);
    vel(2) += sum(f_.col(2) % coeffs);
    coeffs = f_.col(0) % v0 + f_.col(1) % v1 + f_.col(2) % v2;
    coeffs = coeffs / tmp;
    vel(0) += sum(v0 % coeffs);
    vel(1) += sum(v1 % coeffs);
    vel(2) += sum(v2 % coeffs);
    if ((i == num_grids_ - 1 || i == num_grids_ - 2) && fixed_base_) {
      vel = rowvec{0, 0, 0};
    }
    pos_next_.row(i) = pos_.row(i) + vel * dt_;
  }
}

// find forces for the next iteration in parallel
void MultiFlagella::find_next_force_parallel() {
  f_.zeros();
  boost::thread_group find_force_threads;
  int step = num_grids_ / n_threads_;
  int end;
  for (int i = 0; i < num_grids_; i += step) {
    end = (i + step > num_grids_) ? num_grids_ : i + step;
    find_force_threads.create_thread(boost::bind(find_force_in_range, this, i, end));
  }
  find_force_threads.join_all();
}

// find the force in a range of grids
void MultiFlagella::find_force_in_range(int start, int end) {
  rowvec v;
  rowvec pos;
  for (int i = start; i < end; i++) {
    pos = pos_.row(i);
    int j = 0;
    double len;
    for (int idx : *(grids_[i]->adj_idx_)) {
      v = pos_.row(idx) - pos;
      len = norm(v);
      f_.row(i) += v * grids_[i]->stiff_const_->at(j) *
                   (len - grids_[i]->rest_len_->at(j)) /
                   grids_[i]->rest_len_->at(j) / len;
      j++;
    }
    if (repulse_flag_ && i < n_ * 3 * 2) {
      int s = n_ * 3, t = n_ * 3 * 2;
      if (i >= n_ * 3) {
        s = 0;
        t = n_ * 3;
      }
      for (int j = s; j < t; j++) {
        v = pos_.row(j) - pos;
        len = norm(v);
        if (len < flagella_min_dist_) {
          f_.row(i) += v * repulseK_ * (len - flagella_min_dist_) /
              flagella_min_dist_ / len;
        }
      }
    }
    if ((i == num_grids_ - 1 || i == num_grids_ - 2) && fixed_base_) {
      f_.row(i) = rowvec{0, 0, 0};
    }
  }
}

// torque direction is along the
void MultiFlagella::update_torque() {
  int offset = n_ * 3;
  for (int i = 0; i < n_flagella_; i++) {
    offset = i * n_ * 3;
    posL_.row(i) = mean(pos_.rows(offset, offset + 2));
    vectL_.row(i) = normalise(cross(pos_.row(offset + 1) - pos_.row(offset),
                                    pos_.row(offset + 2) - pos_.row(offset))) * motorT_;
  }
}

// this is not used just for single flagella case
void MultiFlagella::run() { }

// Run SingleFlagella in multiple threads
void MultiFlagella::run_parallel(int n_threads) {
  n_threads_ = n_threads;
  build_grids_double();
  for (int i = 0; i < num_grids_; i++) {
    double *pos = grids_[i]->pos_;
    rowvec v = {pos[0], pos[1], pos[2]};
    pos_.row(i) = v;
//    Utils::show_row(v);
//    cout << grids_[i]->adj_idx_->size() << endl;
  }

  for (int i = 0; i < num_iters_; i++) {
    iter_ = i;
    update_torque();
    find_next_pos_parallel();
    find_next_force_parallel();

    if (iter_ % 1000 == 0 || iter_ == num_iters_ - 1) {
      cout << "Progress: " << (iter_ + 1.0) / num_iters_ * 100 << "%" << endl;
    }
    if (iter_ % 1000 == 0 || iter_ == num_iters_ - 1) {
      emit_result();
    }
  }

  pos_out_.close();
  f_out_.close();
//  Utils::show_mat(f_);
}

// output results to file
void MultiFlagella::emit_result() {
  for (int i = 0; i < num_grids_; i++) {
    pos_out_ << pos_(i, 0) << " " << pos_(i, 1) << " " << pos_(i, 2) << endl;
    f_out_ << f_(i, 0) << " " << f_(i, 1) << " " << f_(i, 2) << endl;
  }
}

}


