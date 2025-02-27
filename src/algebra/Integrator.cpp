// Copyright (c) Stanford University, The Regents of the University of
//               California, and others.
//
// All Rights Reserved.
//
// See Copyright-SimVascular.txt for additional details.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "Integrator.h"

Integrator::Integrator(Model* model, double time_step_size, double rho,
                       double atol, int max_iter) {
  this->model = model;
  alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho);
  alpha_f = 1.0 / (1.0 + rho);
  gamma = 0.5 + alpha_m - alpha_f;
  ydot_init_coeff = 1.0 - 1.0 / gamma;
  //ydot_init_coeff = 1.0 - 0.5 / gamma;

  y_coeff = gamma * time_step_size;
  y_coeff_jacobian = alpha_f * y_coeff;

  size = model->dofhandler.size();
  system = SparseSystem(size);
  this->time_step_size = time_step_size;
  this->atol = atol;
  this->max_iter = max_iter;

  y_af = Eigen::Matrix<double, Eigen::Dynamic, 1>(size);
  ydot_am = Eigen::Matrix<double, Eigen::Dynamic, 1>(size);

  // Make some memory reservations
  system.reserve(model);
}

// Must declare default constructord and dedtructor
// because of Eigen.
Integrator::Integrator() {}
Integrator::~Integrator() {}

void Integrator::clean() {
  // Cannot be in destructor because dynamically allocated pointers will be lost
  // when objects are assigned from temporary objects.
  system.clean();
}

void Integrator::update_params(double time_step_size) {
  this->time_step_size = time_step_size;
  y_coeff = gamma * time_step_size;
  y_coeff_jacobian = alpha_f * y_coeff;
  model->update_constant(system);
  model->update_time(system, 0.0);
}

State Integrator::step(const State& old_state, double time) {
  // Predictor: Constant y, consistent ydot
  State new_state = State::Zero(size);
  new_state.ydot += old_state.ydot * ydot_init_coeff;
  //new_state.y += old_state.y;
  new_state.y += old_state.y + 0.5 * old_state.ydot * time_step_size;

  // Determine new time (evaluate terms at generalized mid-point)
  double new_time = time + alpha_f * time_step_size;

  // Evaluate time-dependent element contributions in system
  model->update_time(system, new_time);

  // Count total number of step calls
  n_iter++;
  std::cout << "Time:" << time << std::endl;
  // Non-linear Newton-Raphson iterations
  for (size_t i = 0; i < max_iter; i++) {
    // Initiator: Evaluate the iterates at the intermediate time levels
    ydot_am.setZero();
    y_af.setZero();
    ydot_am += (old_state.ydot + (new_state.ydot - old_state.ydot) * alpha_m);
    y_af += old_state.y + (new_state.y - old_state.y) * alpha_f;
    //std::cout << "ydot_am: " << ydot_am << std::endl;
    //std::cout << "y_af: " << y_af << std::endl;
    // Update solution-dependent element contribitions
    model->update_solution(system, y_af, ydot_am);

    // Evaluate residual
    system.update_residual(y_af, ydot_am);
    
    std::cout << "Residual: " << system.residual.cwiseAbs().maxCoeff() << std::endl;
    // Check termination criterium
    if (system.residual.cwiseAbs().maxCoeff() < atol) {
      std::cout << "Advancing in time with residual" << system.residual.cwiseAbs().maxCoeff() << std::endl;
      break;
    }

    // Abort if maximum number of non-linear iterations is reached
    // else if (i >= max_iter - 20) {
    //   system.get_cond();
    // }
    else if (i == max_iter - 1) {
      std::cout << "Maximum number of non-linear iterations :" << max_iter << std::endl;
      std::cout << "time step size: " << time_step_size << std::endl;
      std::cout << "alpha_f: " << alpha_f << std::endl;
      std::cout << "atol: " << atol << std::endl;
      throw std::runtime_error(
          "Maximum number of non-linear iterations reached.");
    }

    // Evaluate Jacobian
    system.update_jacobian(alpha_m, y_coeff_jacobian);

    // Solve system for increment in ydot
    system.solve();

    // Line Search
    double current_residual = system.residual.cwiseAbs().maxCoeff();
    double new_residual = current_residual;

    bool use_line_search = false;
    double step_size = 1.0;
    double max_num_ls = 10;
    double num_ls = 0;

    while (new_residual >= current_residual) {
      

      State hyp_state = State::Zero(size);
      hyp_state.y += new_state.y;
      hyp_state.ydot += new_state.ydot;

      if (num_ls == 20) {
        step_size = -1;
      }

      hyp_state.ydot += system.dydot * step_size;
      hyp_state.y += system.dydot * step_size * y_coeff;

      ydot_am.setZero();
      y_af.setZero();
      ydot_am += (old_state.ydot + (hyp_state.ydot - old_state.ydot) * alpha_m);
      y_af += old_state.y + (hyp_state.y - old_state.y) * alpha_f;

      // Update solution-dependent element contribitions
      model->update_solution(system, y_af, ydot_am);

      // Evaluate residual
      system.update_residual(y_af, ydot_am);
      new_residual = system.residual.cwiseAbs().maxCoeff();
      if (use_line_search == true){
          step_size *= 0.5;
          num_ls += 1;
      } else {
          break;
      }

      DEBUG_MSG("Step size: " << step_size << " Residual: " << new_residual);
      // Newton's method stuck in local minimum.  Randomize solution and restart
      if (num_ls == max_num_ls * 2) {
        double min_range = 0.0;
        double max_range = 5;
        Eigen::MatrixXd rand = Eigen::MatrixXd::Random(old_state.y.rows(),1);
        rand = (rand + Eigen::MatrixXd::Constant(old_state.y.rows(),1,1.0)) / 2.0;

        State old_state = State::Zero(size);
        old_state.y += rand * (max_range - min_range) + Eigen::MatrixXd::Constant(old_state.y.rows(),1,min_range);
        old_state.ydot += old_state.ydot * 0;

        State new_state = State::Zero(size);
        new_state.ydot += old_state.ydot * ydot_init_coeff;
        new_state.y += old_state.y + 0.5 * old_state.ydot * time_step_size;

        ydot_am.setZero();
        y_af.setZero();
        ydot_am += (old_state.ydot + (new_state.ydot - old_state.ydot) * alpha_m);
        y_af += old_state.y + (new_state.y - old_state.y) * alpha_f;
          
        model->update_solution(system, y_af, ydot_am);
        system.update_residual(y_af, ydot_am);
        system.solve();

        // Reset line search
        current_residual = system.residual.cwiseAbs().maxCoeff();
        new_residual = system.residual.cwiseAbs().maxCoeff();
        step_size = 1.0;
        num_ls = 0;
        DEBUG_MSG("Randomized solution, restarting Newton's method");
      }
    }

    // End of line search

    // Perform post-solve actions on blocks
    model->post_solve(new_state.y);

    // Update the solution
    new_state.ydot += system.dydot * step_size;
    new_state.y += system.dydot * step_size * y_coeff;


    // Count total number of nonlinear iterations
    n_nonlin_iter++;
  }

  return new_state;
}

double Integrator::avg_nonlin_iter() {
  return (double)n_nonlin_iter / (double)n_iter;
}
