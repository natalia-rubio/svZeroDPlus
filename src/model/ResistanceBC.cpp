// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ResistanceBC.h"
#include "../solve/debug.h"

void ResistanceBC::setup_dofs(DOFHandler& dofhandler) {
  Block::setup_dofs_(dofhandler, 1, {});
}

void ResistanceBC::update_constant(SparseSystem& system,
                                   std::vector<double>& parameters) {
  DEBUG_MSG("[ResistanceBC::update_constant] Starting");
  DEBUG_MSG("[ResistanceBC::update_constant] global_eqn_ids size: " << global_eqn_ids.size());
  DEBUG_MSG("[ResistanceBC::update_constant] global_var_ids size: " << global_var_ids.size());
  if (global_eqn_ids.size() > 0) {
    DEBUG_MSG("[ResistanceBC::update_constant] global_eqn_ids[0] = " << global_eqn_ids[0]);
  } else {
    DEBUG_MSG("[ResistanceBC::update_constant] ERROR: global_eqn_ids is empty!");
  }
  if (global_var_ids.size() > 0) {
    DEBUG_MSG("[ResistanceBC::update_constant] global_var_ids[0] = " << global_var_ids[0]);
  } else {
    DEBUG_MSG("[ResistanceBC::update_constant] ERROR: global_var_ids is empty!");
  }
  DEBUG_MSG("[ResistanceBC::update_constant] system.F size: " << system.F.rows() << "x" << system.F.cols());
  DEBUG_MSG("[ResistanceBC::update_constant] About to set F coefficient...");
  if (global_eqn_ids.size() > 0 && global_var_ids.size() > 0) {
    if (global_eqn_ids[0] >= 0 && global_eqn_ids[0] < system.F.rows() &&
        global_var_ids[0] >= 0 && global_var_ids[0] < system.F.cols()) {
      system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
      DEBUG_MSG("[ResistanceBC::update_constant] F coefficient set successfully");
    } else {
      DEBUG_MSG("[ResistanceBC::update_constant] ERROR: Indices out of bounds!");
      DEBUG_MSG("[ResistanceBC::update_constant]   global_eqn_ids[0]=" << global_eqn_ids[0] << ", F.rows()=" << system.F.rows());
      DEBUG_MSG("[ResistanceBC::update_constant]   global_var_ids[0]=" << global_var_ids[0] << ", F.cols()=" << system.F.cols());
    }
  }
  DEBUG_MSG("[ResistanceBC::update_constant] Completed");
}

void ResistanceBC::update_time(SparseSystem& system,
                               std::vector<double>& parameters) {
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) =
      -parameters[global_param_ids[0]];
  system.C(global_eqn_ids[0]) = -parameters[global_param_ids[1]];
}
