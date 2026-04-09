// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "calibrate.h"

#include <limits>
#include <map>

#include "LevenbergMarquardtOptimizer.h"
#include "SimulationParameters.h"

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
  bool calibrate_capacitance =
      calibration_parameters.value("calibrate_capacitance", false);
  double lambda0 = calibration_parameters.value("initial_damping_factor", 1.0);
  double l2_penalty_R =
      calibration_parameters.value("L2_penalty_R_poiseuille", 0.0);
  double l2_penalty_stenosis =
      calibration_parameters.value("L2_penalty_stenosis_coefficient", 0.0);
  double l2_penalty_L =
      calibration_parameters.value("L2_penalty_L", 0.0);

  // Store initial capacitance values (to restore if not calibrating capacitance)
  std::map<std::string, double> initial_capacitance;
  
  // Print Penalty values
  std::cout << "L2 penalty on R_poiseuille: weight = " << l2_penalty_R << std::endl;
  std::cout << "L2 penalty on stenosis_coefficient: weight = " << l2_penalty_stenosis << std::endl;
  std::cout << "L2 penalty on L: weight = " << l2_penalty_L << std::endl;
  // Print capacitance calibration status
  if (zero_capacitance) {
    std::cout << "Capacitance: Setting all to zero" << std::endl;
  } else if (!calibrate_capacitance) {
    std::cout << "Capacitance: Using BloodVesselFC blocks (fixed capacitance, not in optimization)" << std::endl;
  } else {
    std::cout << "Capacitance: Including in optimization (BloodVessel blocks)" << std::endl;
  }

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
  
  // Determine number of optimization parameters per vessel
  // BloodVessel: R, C, L, [stenosis] -> 3 or 4 params
  // BloodVesselFC: R, L, [stenosis] -> 2 or 3 params (no C)
  int num_vessel_params = num_params;
  if (!calibrate_capacitance && !zero_capacitance) {
    num_vessel_params = num_params - 1;  // Skip capacitance parameter
  }
  
  for (auto const& vessel_config : config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];

    // Read capacitance from config and store if using fixed capacitance
    if (!calibrate_capacitance && !zero_capacitance) {
      double cap_value = vessel_config["zero_d_element_values"]["C"].get<double>();
      model.fixed_capacitance[vessel_name] = cap_value;
      initial_capacitance[vessel_name] = cap_value;
    }

    bool is_connector = (vessel_name.find("connector") != std::string::npos);

    // Create parameter IDs
    std::vector<int> param_ids;
    for (size_t k = 0; k < num_vessel_params; k++) {
      int param_id = param_counter++;
      param_ids.push_back(param_id);
      if (is_connector) {
        fixed_param_ids.push_back(param_id);
      }
    }
    if (is_connector) {
      std::cout << "Fixing vessel params for connector: " << vessel_name << std::endl;
    }
    
    // Choose block type: BloodVesselFC for fixed capacitance, otherwise original
    std::string block_type;
    if (!calibrate_capacitance && !zero_capacitance) {
      block_type = "BloodVesselFC";
    } else {
      block_type = vessel_config["zero_d_element_type"].get<std::string>();
    }
    
    model.add_block(block_type, param_ids, vessel_name);
    vessel_id_map.insert({vessel_config["vessel_id"], vessel_name});
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

  // Create junctions
  int num_junc_params_per_outlet = num_params - 1;  // R, L, stenosis (no C)
  for (auto const& junction_config : config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];
    auto const& outlet_vessels = junction_config["outlet_vessels"];
    int num_outlets = outlet_vessels.size();
    
    // Read junction type from input (default to NORMAL_JUNCTION)
    std::string junction_type = junction_config.value("junction_type", "NORMAL_JUNCTION");

    if (num_outlets == 1) {
      // Single outlet junctions are always NORMAL_JUNCTION
      model.add_block("NORMAL_JUNCTION", {}, junction_name);
    } else if (junction_type == "NORMAL_JUNCTION") {
      // Multi-outlet NORMAL_JUNCTION: no parameters to calibrate
      model.add_block("NORMAL_JUNCTION", {}, junction_name);
    } else {
      // Multi-outlet BloodVesselJunction (or other types): has R, L, stenosis parameters
      int junc_param_start = param_counter;
      std::vector<int> param_ids;
      for (size_t i = 0; i < (num_outlets * num_junc_params_per_outlet); i++)
        param_ids.push_back(param_counter++);
      model.add_block("BloodVesselJunction", param_ids, junction_name);

      // Fix junction params for outlets connected to non-EL connector vessels.
      // Parameter layout: [R_0..R_{n-1}, L_0..L_{n-1}, S_0..S_{n-1}]
      for (int oi = 0; oi < num_outlets; oi++) {
        int64_t ov_id = outlet_vessels[oi].get<int64_t>();
        std::string ov_name = vessel_id_map[ov_id];
        bool is_connector = (ov_name.find("connector") != std::string::npos);
        bool is_el_connector = (ov_name.find("connectorEL") != std::string::npos);
        if (is_connector && !is_el_connector) {
          fixed_param_ids.push_back(junc_param_start + oi);                      // R
          fixed_param_ids.push_back(junc_param_start + num_outlets + oi);        // L
          if (num_params > 3) {
            fixed_param_ids.push_back(junc_param_start + 2 * num_outlets + oi);  // stenosis
          }
          std::cout << "Fixing junction " << junction_name
                    << " outlet params for non-EL connector: " << ov_name << std::endl;
        }
      }
    }

    // Check for connections to inlet and outlet vessels and append to
    // connections list
    for (auto vessel_id : junction_config["inlet_vessels"]) {
      connections.push_back({vessel_id_map[vessel_id], junction_name});
    }

    for (auto vessel_id : outlet_vessels) {
      connections.push_back({junction_name, vessel_id_map[vessel_id]});
    }
    DEBUG_MSG("Created junction " << junction_name << " (" << junction_type << ")");
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

  DEBUG_MSG("Number of parameters " << param_counter);

  // Read observations
  DEBUG_MSG("Reading observations");
  int num_obs = 0;
  std::vector<std::vector<double>> y_all;
  std::vector<std::vector<double>> dy_all;
  auto y_values = config["y"];
  auto dy_values = config["dy"];
  
  // Determine num_obs from first available observation
  for (auto& [key, val] : y_values.items()) {
    num_obs = val.get<std::vector<double>>().size();
    break;
  }
  if (num_obs == 0) {
    std::cout << "ERROR: No observations found in 'y'" << std::endl;
    exit(1);
  }
  
  // Initialize observation vectors
  y_all.resize(num_obs);
  dy_all.resize(num_obs);
  
  const double NaN_VALUE = std::numeric_limits<double>::quiet_NaN();
  int missing_count = 0;
  
  for (size_t i = 0; i < model.dofhandler.get_num_variables(); i++) {
    std::string var_name = model.dofhandler.variables[i];
    DEBUG_MSG("Reading observations for variable " << var_name);
    
    if (y_values.contains(var_name) && dy_values.contains(var_name)) {
      auto y_array = y_values[var_name].get<std::vector<double>>();
      auto dy_array = dy_values[var_name].get<std::vector<double>>();
      
      if (y_array.size() != num_obs || dy_array.size() != num_obs) {
        std::cout << "WARNING: Observation size mismatch for '" << var_name 
                  << "'. Expected " << num_obs << " points, got y:" 
                  << y_array.size() << " dy:" << dy_array.size() << std::endl;
      }
      
      for (size_t j = 0; j < num_obs; j++) {
        if (j < y_array.size()) {
          y_all[j].push_back(y_array[j]);
        } else {
          y_all[j].push_back(NaN_VALUE);
        }
        if (j < dy_array.size()) {
          dy_all[j].push_back(dy_array[j]);
        } else {
          dy_all[j].push_back(NaN_VALUE);
        }
      }
    } else {
      // Missing observation - fill with NaN
      missing_count++;
      std::cout << "WARNING: Missing observation for '" << var_name 
                << "'. Using NaN placeholder." << std::endl;
      for (size_t j = 0; j < num_obs; j++) {
        y_all[j].push_back(NaN_VALUE);
        dy_all[j].push_back(NaN_VALUE);
      }
    }
  }
  
  if (missing_count > 0) {
    std::cout << "INFO: " << missing_count << " variable(s) missing observations. "
              << "Residual computation will skip these variables." << std::endl;
  }
  DEBUG_MSG("Number of observations: " << num_obs);

  // Build param index lists for L2 penalty on R_poiseuille, L, and stenosis_coefficient
  std::vector<int> r_poiseuille_param_ids;
  std::vector<int> l_param_ids;
  std::vector<int> stenosis_param_ids;
  if (l2_penalty_R != 0.0 || l2_penalty_L != 0.0 || l2_penalty_stenosis != 0.0) {
    for (auto& vessel_config : output_config["vessels"]) {
      std::string vessel_name = vessel_config["vessel_name"];
      auto block = model.get_block(vessel_name);
      if (!calibrate_capacitance && !zero_capacitance) {
        // BloodVesselFC: R(0), L(1), stenosis(2)
        r_poiseuille_param_ids.push_back(block->global_param_ids[0]);
        l_param_ids.push_back(block->global_param_ids[1]);
        if (calibrate_stenosis) {
          stenosis_param_ids.push_back(block->global_param_ids[2]);
        }
      } else {
        // BloodVessel: R(0), C(1), L(2), stenosis(3)
        r_poiseuille_param_ids.push_back(block->global_param_ids[0]);
        l_param_ids.push_back(block->global_param_ids[2]);
        if (num_params > 3) {
          stenosis_param_ids.push_back(block->global_param_ids[3]);
        }
      }
    }
    for (auto& junction_config : output_config["junctions"]) {
      std::string junction_name = junction_config["junction_name"];
      auto block = model.get_block(junction_name);
      int num_outlets = block->outlet_nodes.size();
      if (num_outlets < 2 || block->global_param_ids.empty()) {
        continue;
      }
      for (int i = 0; i < num_outlets; i++) {
        r_poiseuille_param_ids.push_back(block->global_param_ids[i]);
        l_param_ids.push_back(block->global_param_ids[i + num_outlets]);
        if (num_params > 3) {
          stenosis_param_ids.push_back(
              block->global_param_ids[i + 2 * num_outlets]);
        }
      }
    }
    if (l2_penalty_R != 0.0) {
      std::cout << "L2 penalty on R_poiseuille: weight = " << l2_penalty_R
                << " (" << r_poiseuille_param_ids.size() << " params)"
                << std::endl;
    }
    if (l2_penalty_L != 0.0) {
      std::cout << "L2 penalty on L: weight = " << l2_penalty_L << " ("
                << l_param_ids.size() << " params)" << std::endl;
    }
    if (l2_penalty_stenosis != 0.0) {
      std::cout << "L2 penalty on stenosis_coefficient: weight = "
                << l2_penalty_stenosis << " (" << stenosis_param_ids.size()
                << " params)" << std::endl;
    }
  }

  // Setup start parameter vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha =
      Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(param_counter);
  DEBUG_MSG("Reading initial alpha");
  for (auto& vessel_config : output_config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];
    DEBUG_MSG("Reading initial alpha for " << vessel_name);
    auto block = model.get_block(vessel_name);
    
    if (!calibrate_capacitance && !zero_capacitance) {
      // BloodVesselFC: params are R(0), L(1), stenosis(2)
      alpha[block->global_param_ids[0]] =
          vessel_config["zero_d_element_values"].value("R_poiseuille", 0.0);
      alpha[block->global_param_ids[1]] =
          vessel_config["zero_d_element_values"].value("L", 0.0);
      if (calibrate_stenosis) {
        alpha[block->global_param_ids[2]] =
            vessel_config["zero_d_element_values"].value("stenosis_coefficient", 0.0);
      }
    } else {
      // BloodVessel: params are R(0), C(1), L(2), stenosis(3)
      alpha[block->global_param_ids[0]] =
          vessel_config["zero_d_element_values"].value("R_poiseuille", 0.0);
      double c_init = vessel_config["zero_d_element_values"].value("C", 0.0);
      alpha[block->global_param_ids[1]] = c_init;
      // Store initial capacitance for later restoration if not calibrating
      initial_capacitance[vessel_name] = c_init;
      alpha[block->global_param_ids[2]] =
          vessel_config["zero_d_element_values"].value("L", 0.0);
      if (num_params > 3) {
        alpha[block->global_param_ids[3]] =
            vessel_config["zero_d_element_values"].value("stenosis_coefficient", 0.0);
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
    
    // Skip junctions without parameters (e.g., NORMAL_JUNCTION)
    if (block->global_param_ids.empty()) {
      continue;
    }

    // Initialize junction parameters to zero
    for (size_t i = 0; i < num_outlets; i++) {
      alpha[block->global_param_ids[i]] = 0.0;
      alpha[block->global_param_ids[i + num_outlets]] = 0.0;
      if (num_params > 3) {
        alpha[block->global_param_ids[i + 2 * num_outlets]] = 0.0;
      }
    }
    
    // Read initial values from BloodVesselJunction if available
    if (junction_config.contains("junction_type") && 
        junction_config["junction_type"] == "BloodVesselJunction" &&
        junction_config.contains("junction_values")) {
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
  }

  // Run optimization
  DEBUG_MSG("Start optimization");
  std::cout << "Fixed " << fixed_param_ids.size() << " parameter(s) out of "
            << param_counter << " total" << std::endl;
  auto lm_alg = LevenbergMarquardtOptimizer(
      &model, num_obs, param_counter, lambda0, gradient_tol, increment_tol,
      max_iter, fixed_param_ids, l2_penalty_R, l2_penalty_L,
      l2_penalty_stenosis, r_poiseuille_param_ids, l_param_ids,
      stenosis_param_ids);

  alpha = lm_alg.run(alpha, y_all, dy_all);

  // Write optimized simulation config file
  for (auto& vessel_config : output_config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];
    auto block = model.get_block(vessel_name);
    
    double r_value, l_value, c_value, stenosis_coeff;
    
    if (!calibrate_capacitance && !zero_capacitance) {
      // BloodVesselFC: params are R(0), L(1), stenosis(2)
      r_value = alpha[block->global_param_ids[0]];
      l_value = alpha[block->global_param_ids[1]];
      c_value = initial_capacitance[vessel_name];
      stenosis_coeff = calibrate_stenosis ? alpha[block->global_param_ids[2]] : 0.0;
      
      // Print l_value for each vessel
      std::cout << "l_value for " << vessel_name << " is " << l_value << std::endl;
    } else {
      // BloodVessel: params are R(0), C(1), L(2), stenosis(3)
      r_value = alpha[block->global_param_ids[0]];
      l_value = alpha[block->global_param_ids[2]];
      stenosis_coeff = (num_params > 3) ? alpha[block->global_param_ids[3]] : 0.0;
      
      // Determine capacitance value
      if (zero_capacitance) {
        c_value = 0.0;
      } else {
        c_value = alpha[block->global_param_ids[1]];
        // Warn if negative
        if (c_value < 0.0) {
          std::cout << "WARNING: Optimized C was " << c_value 
                    << " for vessel " << vessel_name 
                    << ". Using absolute value." << std::endl;
          c_value = std::abs(c_value);
        }
      }
      
      // Warn if L is negative
      // if (l_value < 0.0) {
      //   std::cout << "WARNING: L was " << l_value 
      //             << " and is being set to zero for vessel " << vessel_name << std::endl;
      //   l_value = 0.0;
      // }
    }
    
    vessel_config["zero_d_element_values"] = {
        {"R_poiseuille", r_value},
        {"C", c_value},
        {"L", l_value},
        {"stenosis_coefficient", stenosis_coeff}};
  }
  for (auto& junction_config : output_config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];
    auto block = model.get_block(junction_name);
    int num_outlets = block->outlet_nodes.size();
    
    // Get junction type from input (preserve it in output)
    std::string junction_type = junction_config.value("junction_type", "NORMAL_JUNCTION");

    if (num_outlets < 2) {
      // Single outlet junctions stay as-is
      continue;
    }
    
    // Check if this junction has parameters to calibrate
    // NORMAL_JUNCTION has no parameters, BloodVesselJunction has parameters
    if (block->global_param_ids.empty()) {
      // No parameters - this is a NORMAL_JUNCTION, preserve it
      // junction_type is already set from input, no junction_values needed
      continue;
    }

    // BloodVesselJunction with calibrated parameters
    std::vector<double> r_values;
    for (size_t i = 0; i < num_outlets; i++) {
      r_values.push_back(alpha[block->global_param_ids[i]]);
    }
    std::vector<double> l_values;
    for (size_t i = 0; i < num_outlets; i++) {
      l_values.push_back(alpha[block->global_param_ids[i + num_outlets]]);
    }

    std::vector<double> ste_values;

    if (num_params > 3) {
      for (size_t i = 0; i < num_outlets; i++) {
        ste_values.push_back(
            alpha[block->global_param_ids[i + 2 * num_outlets]]);
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
