// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "calibrate.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>

#include "LevenbergMarquardtOptimizer.h"
#include "SimulationParameters.h"

namespace {

// Match naming used for Windkessel connector segments (exclude elastance lumped
// "connectorEL" paths).
bool vessel_is_non_el_connector(const std::string& vessel_name) {
  return vessel_name.find("connector") != std::string::npos &&
         vessel_name.find("connectorEL") == std::string::npos;
}

/** svZeroD JSON may store ids as float (e.g. 2.0); nlohmann get<int64_t> throws type_error. */
std::int64_t json_to_int64(const nlohmann::json& j) {
  if (j.is_number_integer()) {
    return j.get<std::int64_t>();
  }
  if (j.is_number_unsigned()) {
    return static_cast<std::int64_t>(j.get<std::uint64_t>());
  }
  if (j.is_number_float()) {
    return static_cast<std::int64_t>(std::llround(j.get<double>()));
  }
  throw std::runtime_error(std::string("json_to_int64: expected JSON number, got ") +
                           j.type_name());
}

int junction_outlet_count(const nlohmann::json& junction_config) {
  if (junction_config.contains("outlet_blocks") &&
      junction_config["outlet_blocks"].is_array() &&
      !junction_config["outlet_blocks"].empty()) {
    return static_cast<int>(junction_config["outlet_blocks"].size());
  }
  const auto& ov =
      junction_config.value("outlet_vessels", nlohmann::json::array());
  return static_cast<int>(ov.size());
}

/** Vessel or downstream junction name for outlet index ``oi`` (block or id topology). */
std::string outlet_neighbor_name(
    const nlohmann::json& junction_config, int oi,
    const std::map<std::int64_t, std::string>& vessel_id_map) {
  if (junction_config.contains("outlet_blocks") &&
      junction_config["outlet_blocks"].is_array() &&
      static_cast<int>(junction_config["outlet_blocks"].size()) > oi) {
    return junction_config["outlet_blocks"][oi].get<std::string>();
  }
  const auto& ov = junction_config.at("outlet_vessels");
  const std::int64_t ov_id = json_to_int64(ov.at(oi));
  return vessel_id_map.at(ov_id);
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
  std::cout << "[calibrate] freeze_connector_segments="
             << (freeze_connector_segments ? "true" : "false") << std::endl;

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
    vessel_id_map.insert({json_to_int64(vessel_config["vessel_id"]), vessel_name});
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

  for (auto const& vessel_config : config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];
    if (vessel_name.find("connector") == std::string::npos) {
      continue;
    }
    if (!freeze_connector_segments) {
      std::cout << "[calibrate] connector segment \"" << vessel_name
                << "\": NOT FROZEN (freeze_connector_segments is false)"
                << std::endl;
    } else if (vessel_name.find("connectorEL") != std::string::npos) {
      std::cout << "[calibrate] connector segment \"" << vessel_name
                << "\": NOT FROZEN (connectorEL is excluded from freeze)"
                << std::endl;
    } else {
      std::cout << "[calibrate] connector segment \"" << vessel_name
                << "\": FROZEN (R,L,C,stenosis held at 0; junction outlet "
                   "params fixed when applicable)"
                << std::endl;
    }
  }

  // Create junctions
  for (auto const& junction_config : config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];
    const int num_outlets = junction_outlet_count(junction_config);

    if (num_outlets == 1) {
      model.add_block("NORMAL_JUNCTION", {}, junction_name);

    } else {
      const int junc_param_start = param_counter;
      std::vector<int> param_ids;
      for (size_t i = 0; i < (num_outlets * (num_params - 1)); i++)
        param_ids.push_back(param_counter++);
      model.add_block("BloodVesselJunction", param_ids, junction_name);

      for (int oi = 0; oi < num_outlets; oi++) {
        const std::string ov_name =
            outlet_neighbor_name(junction_config, oi, vessel_id_map);
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

    // Check for connections to inlet and outlet vessels and append to
    // connections list (supports inlet_blocks / outlet_blocks for cascaded bifurcations).
    // BloodVesselJunction supports exactly one inlet; use only the first inlet_blocks entry
    // if the JSON lists more than one (avoids duplicate edges / "multiple inlets" at finalize).
    if (junction_config.contains("inlet_blocks") &&
        junction_config["inlet_blocks"].is_array() &&
        !junction_config["inlet_blocks"].empty()) {
      const std::string up =
          junction_config["inlet_blocks"][0].get<std::string>();
      connections.push_back({up, junction_name});
    } else {
      for (auto const& vessel_id :
           junction_config.value("inlet_vessels", nlohmann::json::array())) {
        connections.push_back(
            {vessel_id_map.at(json_to_int64(vessel_id)), junction_name});
      }
    }

    if (junction_config.contains("outlet_blocks") &&
        junction_config["outlet_blocks"].is_array() &&
        !junction_config["outlet_blocks"].empty()) {
      for (auto const& oblock : junction_config["outlet_blocks"]) {
        const std::string dn = oblock.get<std::string>();
        connections.push_back({junction_name, dn});
      }
    } else {
      for (auto const& vessel_id :
           junction_config.value("outlet_vessels", nlohmann::json::array())) {
        connections.push_back(
            {junction_name, vessel_id_map.at(json_to_int64(vessel_id))});
      }
    }
    DEBUG_MSG("Created junction " << junction_name);
  }

  // Drop duplicate directed edges (same upstream -> same downstream) from bad JSON merges.
  {
    std::set<std::pair<std::string, std::string>> conn_seen;
    decltype(connections) conn_deduped;
    conn_deduped.reserve(connections.size());
    for (const auto& c : connections) {
      const auto key =
          std::make_pair(std::get<0>(c), std::get<1>(c));
      if (conn_seen.insert(key).second) {
        conn_deduped.push_back(c);
      }
    }
    connections.swap(conn_deduped);
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
  std::cout << "[calibrate] total optimization parameters: " << param_counter
            << ", fixed parameter indices: " << fixed_param_ids.size()
            << std::endl;

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
      for (int oi = 0; oi < num_outlets; oi++) {
        const std::string ov_name =
            outlet_neighbor_name(junction_config, oi, vessel_id_map);
        if (vessel_is_non_el_connector(ov_name)) {
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

    std::vector<double> r_values;
    for (size_t i = 0; i < static_cast<size_t>(num_outlets); i++) {
      double rv = alpha[block->global_param_ids[i]];
      const std::string ov_name =
          outlet_neighbor_name(junction_config, static_cast<int>(i), vessel_id_map);
      if (freeze_connector_segments &&
          vessel_is_non_el_connector(ov_name)) {
        rv = 0.0;
      }
      r_values.push_back(rv);
    }
    std::vector<double> l_values;
    for (size_t i = 0; i < static_cast<size_t>(num_outlets); i++) {
      double lv =
          std::max(alpha[block->global_param_ids[i + num_outlets]], 0.0);
      const std::string ov_name =
          outlet_neighbor_name(junction_config, static_cast<int>(i), vessel_id_map);
      if (freeze_connector_segments &&
          vessel_is_non_el_connector(ov_name)) {
        lv = 0.0;
      }
      l_values.push_back(lv);
    }

    std::vector<double> ste_values;

    if (num_params > 3) {
      for (size_t i = 0; i < static_cast<size_t>(num_outlets); i++) {
        double sv =
            alpha[block->global_param_ids[i + 2 * num_outlets]];
        const std::string ov_name =
            outlet_neighbor_name(junction_config, static_cast<int>(i), vessel_id_map);
        if (freeze_connector_segments &&
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
