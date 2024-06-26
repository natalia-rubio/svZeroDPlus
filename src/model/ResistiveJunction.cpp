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

#include "ResistiveJunction.h"

void ResistiveJunction::setup_dofs(DOFHandler &dofhandler) {
  // Set number of equations of a junction block based on number of
  // inlets/outlets. Must be set before calling parent constructor
  num_inlets = inlet_nodes.size();
  num_outlets = outlet_nodes.size();
  Block::setup_dofs_(dofhandler, num_inlets + num_outlets + 1, {"pressure_c"});
  num_triplets.F = (num_inlets + num_outlets) * 4;
}

void ResistiveJunction::update_constant(SparseSystem &system,
                                        std::vector<double> &parameters) {
  for (size_t i = 0; i < num_inlets; i++) {
    system.F.coeffRef(global_eqn_ids[i], global_var_ids[i * 2]) = 1.0;
    system.F.coeffRef(global_eqn_ids[i], global_var_ids[i * 2 + 1]) =
        -parameters[global_param_ids[i]];
    system.F.coeffRef(global_eqn_ids[i], global_var_ids.back()) = -1.0;
  }

  for (size_t i = num_inlets; i < num_inlets + num_outlets; i++) {
    system.F.coeffRef(global_eqn_ids[i], global_var_ids[i * 2]) = -1.0;
    system.F.coeffRef(global_eqn_ids[i], global_var_ids[i * 2 + 1]) =
        -parameters[global_param_ids[i]];
    system.F.coeffRef(global_eqn_ids[i], global_var_ids.back()) = 1.0;
  }

  // Mass conservation
  for (size_t i = 1; i < num_inlets * 2; i = i + 2) {
    system.F.coeffRef(global_eqn_ids[num_inlets + num_outlets],
                      global_var_ids[i]) = 1.0;
  }

  for (size_t i = (num_inlets * 2) + 1; i < (num_inlets + num_outlets) * 2;
       i = i + 2) {
    system.F.coeffRef(global_eqn_ids[num_inlets + num_outlets],
                      global_var_ids[i]) = -1.0;
  }
}
