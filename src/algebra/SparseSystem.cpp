// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "SparseSystem.h"

#include "Model.h"
#include "../solve/debug.h"

SparseSystem::SparseSystem() {}

SparseSystem::SparseSystem(int n) {
  DEBUG_MSG("[SparseSystem] Constructor called with n=" << n);
  DEBUG_MSG("[SparseSystem] Creating sparse matrices...");
  F = Eigen::SparseMatrix<double>(n, n);
  E = Eigen::SparseMatrix<double>(n, n);
  dC_dy = Eigen::SparseMatrix<double>(n, n);
  dC_dydot = Eigen::SparseMatrix<double>(n, n);
  DEBUG_MSG("[SparseSystem] Sparse matrices created");
  C = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n);
  DEBUG_MSG("[SparseSystem] C vector created");

  jacobian = Eigen::SparseMatrix<double>(n, n);
  residual = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n);
  dydot = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n);
  DEBUG_MSG("[SparseSystem] Constructor completed");
}

SparseSystem::~SparseSystem() {}

void SparseSystem::clean() {
  // Cannot be in destructor because dynamically allocated pointers will be lost
  // when objects are assigned from temporary objects.
  // delete solver;
}

void SparseSystem::reserve(Model* model) {
  DEBUG_MSG("[SparseSystem::reserve] Starting reserve");
  DEBUG_MSG("[SparseSystem::reserve] Getting num_triplets from model...");
  auto num_triplets = model->get_num_triplets();
  DEBUG_MSG("[SparseSystem::reserve] num_triplets: F=" << num_triplets.F << ", E=" << num_triplets.E << ", D=" << num_triplets.D);
  DEBUG_MSG("[SparseSystem::reserve] Reserving memory for sparse matrices...");
  F.reserve(num_triplets.F);
  E.reserve(num_triplets.E);
  dC_dy.reserve(num_triplets.D);
  dC_dydot.reserve(num_triplets.D);
  DEBUG_MSG("[SparseSystem::reserve] Memory reserved");

  DEBUG_MSG("[SparseSystem::reserve] Calling model->update_constant...");
  model->update_constant(*this);
  DEBUG_MSG("[SparseSystem::reserve] model->update_constant completed");
  DEBUG_MSG("[SparseSystem::reserve] Calling model->update_time...");
  model->update_time(*this, 0.0);
  DEBUG_MSG("[SparseSystem::reserve] model->update_time completed");

  DEBUG_MSG("[SparseSystem::reserve] Creating dummy vectors...");
  Eigen::Matrix<double, Eigen::Dynamic, 1> dummy_y =
      Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(residual.size());

  Eigen::Matrix<double, Eigen::Dynamic, 1> dummy_dy =
      Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(residual.size());
  DEBUG_MSG("[SparseSystem::reserve] Dummy vectors created");

  DEBUG_MSG("[SparseSystem::reserve] Calling model->update_solution...");
  model->update_solution(*this, dummy_y, dummy_dy);
  DEBUG_MSG("[SparseSystem::reserve] model->update_solution completed");

  F.makeCompressed();
  E.makeCompressed();
  dC_dy.makeCompressed();
  dC_dydot.makeCompressed();
  jacobian.reserve(num_triplets.F + num_triplets.E);  // Just an estimate
  update_jacobian(1.0, 1.0);  // Update it once to have sparsity pattern
  jacobian.makeCompressed();
  solver->analyzePattern(jacobian);  // Let solver analyze pattern
}

void SparseSystem::update_residual(
    Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    Eigen::Matrix<double, Eigen::Dynamic, 1>& ydot) {
  residual.setZero();
  residual -= C;
  residual.noalias() -= E * ydot;
  residual.noalias() -= F * y;
}

void SparseSystem::update_jacobian(double time_coeff_ydot,
                                   double time_coeff_y) {
  jacobian.setZero();
  jacobian += (E + dC_dydot) * time_coeff_ydot;
  jacobian += (F + dC_dy) * time_coeff_y;
}

void SparseSystem::solve() {
  solver->factorize(jacobian);
  if (solver->info() != Eigen::Success) {
    throw std::runtime_error(
        "System is singular. Check your model (connections, boundary "
        "conditions, parameters).");
  }
  dydot.setZero();
  dydot += solver->solve(residual);
}
