// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "calibrate.h"

#include <algorithm>
#include <iostream>
#include <unordered_set>

#include "LevenbergMarquardtOptimizer.h"
#include "SimulationParameters.h"

namespace {

// Match naming used for Windkessel connector segments (exclude elastance lumped
// "connectorEL" paths).
bool vessel_is_non_el_connector(const std::string& vessel_name) {
  return vessel_name.find("connector") != std::string::npos &&
         vessel_name.find("connectorEL") == std::string::npos;
}

}  // namespace

nlohmann::json calibrate(const nlohmann::json& config) {
  auto output_config = nlohmann::json(config);

  // Read calibration parameters
  DEBUG_MSG("Parse calibration parameters");
  auto const& calibration_parameters = config["calibration_parameters"];
  double gradient_tol =
      calibration_parameters.value("tolerance_gradient", 1e-5);
  double increment_tol =
      calibration_parameters.value("tolerance_increment", 1e-10);
  int max_iter = calibration_parameters.value("maximum_iterations", 100);
  bool calibrate_stenosis =
      calibration_parameters.value("calibrate_stenosis_coefficient", true);
  bool zero_capacitance =
      calibration_parameters.value("set_capacitance_to_zero", false);
  bool freeze_connector_segments =
      calibration_parameters.value("freeze_connector_segments", true);
  double lambda0 = calibration_parameters.value("initial_damping_factor", 1.0);

  int num_params = 3;
  if (calibrate_stenosis) {
    num_params = 4;
  }

  // Setup model
  auto model = Model();
  std::vector<std::tuple<std::string, std::string>> connections;
  std::vector<std::tuple<std::string, std::string>> inlet_connections;
  std::vector<std::tuple<std::string, std::string>> outlet_connections;

  // Create vessels
  DEBUG_MSG("Load vessels");
  std::map<std::int64_t, std::string> vessel_id_map;
  int param_counter = 0;
  std::vector<int> fixed_param_ids;
  for (auto const& vessel_config : config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];

    std::string block_type =
        vessel_config["zero_d_element_type"].get<std::string>();
    const int vessel_param_start = param_counter;
    std::vector<int> param_ids;
    if (block_type == "BloodVesselFC") {
      int n_fc = calibrate_stenosis ? 3 : 2;
      for (int k = 0; k < n_fc; k++) {
        param_ids.push_back(param_counter++);
      }
    } else {
      for (size_t k = 0; k < num_params; k++) {
        param_ids.push_back(param_counter++);
      }
    }
    if (freeze_connector_segments &&
        vessel_is_non_el_connector(vessel_name)) {
      if (block_type == "BloodVesselFC") {
        fixed_param_ids.push_back(vessel_param_start);      // R
        fixed_param_ids.push_back(vessel_param_start + 1);  // L
        if (calibrate_stenosis) {
          fixed_param_ids.push_back(vessel_param_start + 2);  // stenosis
        }
      } else {
        fixed_param_ids.push_back(vessel_param_start);      // R
        fixed_param_ids.push_back(vessel_param_start + 1);  // C
        fixed_param_ids.push_back(vessel_param_start + 2);  // L
        if (num_params > 3) {
          fixed_param_ids.push_back(vessel_param_start + 3);  // stenosis
        }
      }
      DEBUG_MSG("Freeze connector vessel " << vessel_name << " (" << block_type
                                            << ")");
    }
    model.add_block(block_type, param_ids, vessel_name);
    vessel_id_map.insert({vessel_config["vessel_id"], vessel_name});
    if (block_type == "BloodVesselFC") {
      double c_val = 0.0;
      if (freeze_connector_segments &&
          vessel_is_non_el_connector(vessel_name)) {
        c_val = 0.0;
      } else if (!zero_capacitance &&
                 vessel_config["zero_d_element_values"].contains("C")) {
        c_val = vessel_config["zero_d_element_values"]["C"].get<double>();
      }
      model.fixed_capacitance[vessel_name] = c_val;
    }
    DEBUG_MSG("Created vessel " << vessel_name);

    // Read connected boundary conditions
    if (vessel_config.contains("boundary_conditions")) {
      auto const& vessel_bc_config = vessel_config["boundary_conditions"];
      if (vessel_bc_config.contains("inlet")) {
        inlet_connections.push_back({vessel_bc_config["inlet"], vessel_name});
      }
      if (vessel_bc_config.contains("outlet")) {
        outlet_connections.push_back({vessel_name, vessel_bc_config["outlet"]});
      }
    }
  }

  // --- Same J-J inlet skip as SimulationParameters::create_junctions ---
  // If inlet_blocks names another junction, we do not add that edge here; the
  // parent junction's outlet_blocks already created (J_parent, J_child).
  std::unordered_set<std::string> junction_names;
  junction_names.reserve(config["junctions"].size());
  for (auto const& jc : config["junctions"]) {
    junction_names.insert(jc["junction_name"].get<std::string>());
  }

  // Create junctions (legacy vessel ids and/or Phase A block lists)
  for (auto const& junction_config : config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];

    // Block connectivity: string names for inlets/outlets (vessels or junctions).
    const bool use_blocks =
        junction_config.contains("inlet_blocks") &&
        junction_config.contains("outlet_blocks") &&
        junction_config["inlet_blocks"].is_array() &&
        junction_config["outlet_blocks"].is_array();
    // Legacy: integer vessel ids on the junction object.
    const bool use_vessel_ids =
        junction_config.contains("inlet_vessels") &&
        junction_config.contains("outlet_vessels");

    // Outlet count drives NORMAL_JUNCTION vs BloodVesselJunction; prefer ids if
    // both representations are present (matches forward-solver JSON rules).
    int num_outlets = 0;
    if (use_vessel_ids) {
      num_outlets = static_cast<int>(junction_config["outlet_vessels"].size());
    } else if (use_blocks) {
      num_outlets = static_cast<int>(junction_config["outlet_blocks"].size());
    }

    if (num_outlets == 1) {
      model.add_block("NORMAL_JUNCTION", {}, junction_name);

    } else if (num_outlets > 1) {
      const int junc_param_start = param_counter;
      std::vector<int> param_ids;
      for (size_t i = 0; i < (num_outlets * (num_params - 1)); i++)
        param_ids.push_back(param_counter++);
      model.add_block("BloodVesselJunction", param_ids, junction_name);

      // Optional: pin R,L,(stenosis) for legacy *connector* vessel outlets.
      if (use_vessel_ids) {
        for (int oi = 0; oi < num_outlets; oi++) {
          std::int64_t ov_id =
              junction_config["outlet_vessels"][oi].get<std::int64_t>();
          const std::string& ov_name = vessel_id_map[ov_id];
          if (freeze_connector_segments &&
              vessel_is_non_el_connector(ov_name)) {
            fixed_param_ids.push_back(junc_param_start + oi);  // R
            fixed_param_ids.push_back(junc_param_start + num_outlets +
                                      oi);  // L
            if (num_params > 3) {
              fixed_param_ids.push_back(junc_param_start + 2 * num_outlets +
                                        oi);  // stenosis
            }
            DEBUG_MSG("Freeze junction " << junction_name
                                          << " outlet R,L"
                                          << (num_params > 3 ? ",stenosis" : "")
                                          << " at 0 for connector " << ov_name);
          }
        }
      }
      // Same freeze policy using outlet block names (vessel names only match
      // vessel_is_non_el_connector; junction outlet names are ignored here).
      else if (use_blocks) {
        for (int oi = 0; oi < num_outlets; oi++) {
          const std::string ov_name =
              junction_config["outlet_blocks"][oi].get<std::string>();
          if (freeze_connector_segments &&
              vessel_is_non_el_connector(ov_name)) {
            fixed_param_ids.push_back(junc_param_start + oi);  // R
            fixed_param_ids.push_back(junc_param_start + num_outlets +
                                      oi);  // L
            if (num_params > 3) {
              fixed_param_ids.push_back(junc_param_start + 2 * num_outlets +
                                        oi);  // stenosis
            }
            DEBUG_MSG("Freeze junction " << junction_name
                                          << " outlet R,L"
                                          << (num_params > 3 ? ",stenosis" : "")
                                          << " at 0 for connector " << ov_name);
          }
        }
      }
    }

    // Graph edges for this junction (prefer vessel ids when both are present).
    if (use_vessel_ids) {
      for (auto vessel_id : junction_config["inlet_vessels"]) {
        connections.push_back({vessel_id_map[vessel_id], junction_name});
      }
      for (auto vessel_id : junction_config["outlet_vessels"]) {
        connections.push_back({junction_name, vessel_id_map[vessel_id]});
      }
    } else if (use_blocks) {
      // Inlets: vessels (and any non-junction block) only; skip J-J duplicates.
      for (const auto& iblock : junction_config["inlet_blocks"]) {
        const std::string up = iblock.get<std::string>();
        // Skip if inlet_blocks names another junction (parent J_parent, child J_child).
        if (junction_names.count(up) != 0) {
          continue;
        }
        connections.push_back({up, junction_name});
      }
      // Outlets: always emit (this junction, downstream block).
      for (const auto& oblock : junction_config["outlet_blocks"]) {
        const std::string dn = oblock.get<std::string>();
        connections.push_back({junction_name, dn});
      }
    }
    DEBUG_MSG("Created junction " << junction_name);
  }

  // Create Connections
  DEBUG_MSG("Created connection");
  for (auto& connection : connections) {
    auto ele1 = model.get_block(std::get<0>(connection));
    auto ele2 = model.get_block(std::get<1>(connection));
    model.add_node({ele1}, {ele2}, ele1->get_name() + ":" + ele2->get_name());
  }
  for (auto& connection : inlet_connections) {
    auto ele = model.get_block(std::get<1>(connection));
    model.add_node({}, {ele}, std::get<0>(connection) + ":" + ele->get_name());
  }
  for (auto& connection : outlet_connections) {
    auto ele = model.get_block(std::get<0>(connection));
    model.add_node({ele}, {}, ele->get_name() + ":" + std::get<1>(connection));
  }

  // Finalize model
  model.finalize();

  std::sort(fixed_param_ids.begin(), fixed_param_ids.end());
  fixed_param_ids.erase(std::unique(fixed_param_ids.begin(), fixed_param_ids.end()),
                        fixed_param_ids.end());
  DEBUG_MSG("Number of parameters " << param_counter << ", fixed indices "
                                    << fixed_param_ids.size());

  // Read observations
  DEBUG_MSG("Reading observations");
  int num_obs = 0;
  std::vector<std::vector<double>> y_all;
  std::vector<std::vector<double>> dy_all;
  auto y_values = config["y"];
  auto dy_values = config["dy"];
  for (size_t i = 0; i < model.dofhandler.get_num_variables(); i++) {
    std::string var_name = model.dofhandler.variables[i];
    DEBUG_MSG("Reading observations for variable " << var_name);
    if (!y_values.contains(var_name)) {
      std::cout << "ERROR: Missing y observation for '" << var_name << "'"
                << std::endl;
      exit(1);
    }
    if (!dy_values.contains(var_name)) {
      std::cout << "ERROR: Missing dy observation for '" << var_name << "'"
                << std::endl;
      exit(1);
    }
    auto y_array = y_values[var_name].get<std::vector<double>>();
    auto dy_array = dy_values[var_name].get<std::vector<double>>();
    num_obs = y_array.size();
    if (i == 0) {
      y_all.resize(num_obs);
      dy_all.resize(num_obs);
    }
    for (size_t j = 0; j < num_obs; j++) {
      y_all[j].push_back(y_array[j]);
      dy_all[j].push_back(dy_array[j]);
    }
  }
  DEBUG_MSG("Number of observations: " << num_obs);

  // Setup start parameter vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha =
      Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(param_counter);
  DEBUG_MSG("Reading initial alpha");
  for (auto& vessel_config : output_config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];
    DEBUG_MSG("Reading initial alpha for " << vessel_name);
    auto block = model.get_block(vessel_name);
    std::string vtype = vessel_config["zero_d_element_type"].get<std::string>();
    if (vtype == "BloodVesselFC") {
      alpha[block->global_param_ids[0]] =
          vessel_config["zero_d_element_values"].value("R_poiseuille", 0.0);
      alpha[block->global_param_ids[1]] =
          vessel_config["zero_d_element_values"].value("L", 0.0);
      if (calibrate_stenosis) {
        alpha[block->global_param_ids[2]] =
            vessel_config["zero_d_element_values"].value("stenosis_coefficient",
                                                         0.0);
      }
    } else {
      alpha[block->global_param_ids[0]] =
          vessel_config["zero_d_element_values"].value("R_poiseuille", 0.0);
      alpha[block->global_param_ids[1]] =
          vessel_config["zero_d_element_values"].value("C", 0.0);
      alpha[block->global_param_ids[2]] =
          vessel_config["zero_d_element_values"].value("L", 0.0);
      if (num_params > 3) {
        alpha[block->global_param_ids[3]] =
            vessel_config["zero_d_element_values"].value("stenosis_coefficient",
                                                         0.0);
      }
    }
    if (freeze_connector_segments &&
        vessel_is_non_el_connector(vessel_name)) {
      if (vtype == "BloodVesselFC") {
        alpha[block->global_param_ids[0]] = 0.0;
        alpha[block->global_param_ids[1]] = 0.0;
        if (calibrate_stenosis) {
          alpha[block->global_param_ids[2]] = 0.0;
        }
      } else {
        alpha[block->global_param_ids[0]] = 0.0;
        alpha[block->global_param_ids[1]] = 0.0;
        alpha[block->global_param_ids[2]] = 0.0;
        if (num_params > 3) {
          alpha[block->global_param_ids[3]] = 0.0;
        }
      }
    }
  }
  for (auto& junction_config : output_config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];
    DEBUG_MSG("Reading initial alpha for " << junction_name);
    auto block = model.get_block(junction_name);
    int num_outlets = block->outlet_nodes.size();

    if (num_outlets < 2) {
      continue;
    }

    for (size_t i = 0; i < num_outlets; i++) {
      alpha[block->global_param_ids[i]] = 0.0;
      alpha[block->global_param_ids[i + num_outlets]] = 0.0;
      if (num_params > 3) {
        alpha[block->global_param_ids[i + 2 * num_outlets]] = 0.0;
      }
    }
    if (junction_config["junction_type"] == "BloodVesselJunction") {
      auto resistance = junction_config["junction_values"]["R_poiseuille"]
                            .get<std::vector<double>>();
      auto inductance =
          junction_config["junction_values"]["L"].get<std::vector<double>>();
      auto stenosis_coeff =
          junction_config["junction_values"]["stenosis_coefficient"]
              .get<std::vector<double>>();
      for (size_t i = 0; i < num_outlets; i++) {
        alpha[block->global_param_ids[i]] = resistance[i];
        alpha[block->global_param_ids[i + num_outlets]] = inductance[i];
        if (num_params > 3) {
          alpha[block->global_param_ids[i + 2 * num_outlets]] =
              stenosis_coeff[i];
        }
      }
    }
    if (freeze_connector_segments && num_outlets >= 2 &&
        !block->global_param_ids.empty()) {
      // Calibration inputs may be legacy id-wired (outlet_vessels) or
      // block-wired (outlet_blocks, including J-J). Resolve a per-outlet name
      // from either representation before applying connector freeze logic.
      // Loop over all outlets and resolve the name from either representation.
      for (int oi = 0; oi < num_outlets; oi++) {
        std::string ov_name;
        // If outlet vessel ids are present, resolve the name from the vessel id map.
        if (junction_config.contains("outlet_vessels") &&
            junction_config["outlet_vessels"].is_array() &&
            static_cast<int>(junction_config["outlet_vessels"].size()) > oi) {
          std::int64_t ov_id =
              junction_config["outlet_vessels"][oi].get<std::int64_t>();
          ov_name = vessel_id_map[ov_id];
        // If outlet block names are present, resolve the name from the outlet block names.
        } else if (junction_config.contains("outlet_blocks") &&
                   junction_config["outlet_blocks"].is_array() &&
                   static_cast<int>(junction_config["outlet_blocks"].size()) >
                       oi) {
          ov_name = junction_config["outlet_blocks"][oi].get<std::string>();
        }
        if (!ov_name.empty() && vessel_is_non_el_connector(ov_name)) {
          alpha[block->global_param_ids[oi]] = 0.0;
          alpha[block->global_param_ids[oi + num_outlets]] = 0.0;
          if (num_params > 3) {
            alpha[block->global_param_ids[oi + 2 * num_outlets]] = 0.0;
          }
        }
      }
    }
  }

  // Run optimization
  const int num_locations =
      static_cast<int>(model.dofhandler.get_num_variables());
  const int num_eq_model = model.dofhandler.get_num_equations();
  const long long num_residual_rows =
      static_cast<long long>(num_obs) * static_cast<long long>(num_eq_model);
  std::cout << "Calibration dimensions: " << num_locations
            << " observed state variables (y keys / DOFs), "
            << num_obs << " time step(s), " << num_eq_model
            << " governing equation(s) per timestep" << std::endl;
  std::cout << "  Stacked residual length: " << num_residual_rows
            << " (= time steps × equations); "
            << "parameters: " << param_counter << std::endl;
  std::cout << "  Jacobian (sparse, least-squares): " << num_residual_rows
            << " x " << param_counter << std::endl;
  std::cout << "  Normal-system matrix (dense LLT each iteration): "
            << param_counter << " x " << param_counter << std::endl;
  DEBUG_MSG("Start optimization");
  auto lm_alg = LevenbergMarquardtOptimizer(
      &model, num_obs, param_counter, lambda0, gradient_tol, increment_tol,
      max_iter, fixed_param_ids);

  alpha = lm_alg.run(alpha, y_all, dy_all);

  // Write optimized simulation config file
  for (auto& vessel_config : output_config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];
    auto block = model.get_block(vessel_name);
    std::string vtype = vessel_config["zero_d_element_type"].get<std::string>();
    const bool freeze_here = freeze_connector_segments &&
                             vessel_is_non_el_connector(vessel_name);

    if (vtype == "BloodVesselFC") {
      double c_value = 0.0;
      if (!zero_capacitance) {
        c_value = model.fixed_capacitance[vessel_name];
      }
      double stenosis_coeff = 0.0;
      if (calibrate_stenosis) {
        stenosis_coeff = alpha[block->global_param_ids[2]];
      }
      double r_out = alpha[block->global_param_ids[0]];
      double l_out = std::max(alpha[block->global_param_ids[1]], 0.0);
      if (freeze_here) {
        r_out = 0.0;
        l_out = 0.0;
        c_value = 0.0;
        stenosis_coeff = 0.0;
      }
      vessel_config["zero_d_element_values"] = {
          {"R_poiseuille", r_out},
          {"C", std::max(c_value, 0.0)},
          {"L", l_out},
          {"stenosis_coefficient", stenosis_coeff}};
    } else {
      double stenosis_coeff = 0.0;
      if (num_params > 3) {
        stenosis_coeff = alpha[block->global_param_ids[3]];
      }
      double c_value = 0.0;
      if (!zero_capacitance) {
        c_value = alpha[block->global_param_ids[1]];
      }
      double r_out = alpha[block->global_param_ids[0]];
      double l_out = std::max(alpha[block->global_param_ids[2]], 0.0);
      if (freeze_here) {
        r_out = 0.0;
        c_value = 0.0;
        l_out = 0.0;
        stenosis_coeff = 0.0;
      }
      vessel_config["zero_d_element_values"] = {
          {"R_poiseuille", r_out},
          {"C", std::max(c_value, 0.0)},
          {"L", l_out},
          {"stenosis_coefficient", stenosis_coeff}};
    }
  }
  for (auto& junction_config : output_config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];
    auto block = model.get_block(junction_name);
    int num_outlets = block->outlet_nodes.size();

    if (num_outlets < 2) {
      continue;
    }

    // Map optimizer outlet index -> output JSON outlet name for both
    // outlet_vessels (ids) and outlet_blocks (names). This keeps
    // freeze_connector_segments behavior consistent for block connectivity.
    auto outlet_name_at = [&](size_t oi) -> std::string {
      if (junction_config.contains("outlet_vessels") &&
          junction_config["outlet_vessels"].is_array() &&
          junction_config["outlet_vessels"].size() > oi) {
        std::int64_t ov_id =
            junction_config["outlet_vessels"][oi].get<std::int64_t>();
        return vessel_id_map[ov_id];
      }
      if (junction_config.contains("outlet_blocks") &&
          junction_config["outlet_blocks"].is_array() &&
          junction_config["outlet_blocks"].size() > oi) {
        return junction_config["outlet_blocks"][oi].get<std::string>();
      }
      return std::string();
    };
    std::vector<double> r_values;
    for (size_t i = 0; i < num_outlets; i++) {
      double rv = alpha[block->global_param_ids[i]];
      const std::string ov_name = outlet_name_at(i);
      if (freeze_connector_segments &&
          !ov_name.empty() &&
          vessel_is_non_el_connector(ov_name)) {
        rv = 0.0;
      }
      r_values.push_back(rv);
    }
    std::vector<double> l_values;
    for (size_t i = 0; i < num_outlets; i++) {
      double lv =
          std::max(alpha[block->global_param_ids[i + num_outlets]], 0.0);
      const std::string ov_name = outlet_name_at(i);
      if (freeze_connector_segments &&
          !ov_name.empty() &&
          vessel_is_non_el_connector(ov_name)) {
        lv = 0.0;
      }
      l_values.push_back(lv);
    }

    std::vector<double> ste_values;

    if (num_params > 3) {
      for (size_t i = 0; i < num_outlets; i++) {
        double sv =
            alpha[block->global_param_ids[i + 2 * num_outlets]];
        const std::string ov_name = outlet_name_at(i);
        if (freeze_connector_segments &&
            !ov_name.empty() &&
            vessel_is_non_el_connector(ov_name)) {
          sv = 0.0;
        }
        ste_values.push_back(sv);
      }
    } else {
      for (size_t i = 0; i < num_outlets; i++) {
        ste_values.push_back(0.0);
      }
    }

    junction_config["junction_type"] = "BloodVesselJunction";
    junction_config["junction_values"] = {{"R_poiseuille", r_values},
                                          {"L", l_values},
                                          {"stenosis_coefficient", ste_values}};
  }

  output_config.erase("y");
  output_config.erase("dy");
  output_config.erase("calibration_parameters");

  return output_config;
}
