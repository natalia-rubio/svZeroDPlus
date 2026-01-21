# Modifying svZeroDCalibrator for Optional Observations

## Current Behavior

The calibrator currently **requires** observations (`y` and `dy`) for every variable in the model. In `calibrate.cpp` (lines 121-145), it:

1. Iterates through all model variables (`model.dofhandler.variables`)
2. Checks if observations exist for each variable (lines 124-133)
3. **Exits with error** if any observations are missing
4. Builds observation matrices `y_all` and `dy_all` where:
   - Each row = one time point
   - Each column = one variable (in order of `model.dofhandler.variables`)

## Proposed Modification Strategy

To make observations optional, we need to modify several parts:

### 1. **Observation Reading (calibrate.cpp, lines 114-146)**

**Current code:**
```cpp
for (size_t i = 0; i < model.dofhandler.get_num_variables(); i++) {
    std::string var_name = model.dofhandler.variables[i];
    if (!y_values.contains(var_name)) {
        std::cout << "ERROR: Missing y observation for '" << var_name << "'" << std::endl;
        exit(1);  // <-- This needs to change
    }
    // ... read observations
}
```

**Modified approach:**
```cpp
// Track which variables have observations
std::vector<bool> has_observation(model.dofhandler.get_num_variables(), false);
std::vector<size_t> observation_indices;  // Maps variable index -> observation column index

for (size_t i = 0; i < model.dofhandler.get_num_variables(); i++) {
    std::string var_name = model.dofhandler.variables[i];
    
    if (y_values.contains(var_name) && dy_values.contains(var_name)) {
        has_observation[i] = true;
        observation_indices.push_back(i);
        
        auto y_array = y_values[var_name].get<std::vector<double>>();
        auto dy_array = dy_values[var_name].get<std::vector<double>>();
        num_obs = y_array.size();
        
        if (observation_indices.size() == 1) {
            y_all.resize(num_obs);
            dy_all.resize(num_obs);
        }
        
        for (size_t j = 0; j < num_obs; j++) {
            y_all[j].push_back(y_array[j]);
            dy_all[j].push_back(dy_array[j]);
        }
    } else {
        // Missing observation - use NaN as placeholder
        if (observation_indices.size() == 0) {
            // Need to determine num_obs from first available observation
            // Or use a default/fallback
        }
        observation_indices.push_back(SIZE_MAX);  // Sentinel value
    }
}
```

**Alternative approach (simpler but less efficient):**
```cpp
// Use NaN as placeholder for missing observations
const double NaN_VALUE = std::numeric_limits<double>::quiet_NaN();

for (size_t i = 0; i < model.dofhandler.get_num_variables(); i++) {
    std::string var_name = model.dofhandler.variables[i];
    
    if (i == 0) {
        // Determine num_obs from first available observation
        if (y_values.contains(var_name)) {
            num_obs = y_values[var_name].get<std::vector<double>>().size();
        } else {
            // Need to find first available observation
            for (auto& [key, val] : y_values.items()) {
                num_obs = val.get<std::vector<double>>().size();
                break;
            }
        }
        y_all.resize(num_obs);
        dy_all.resize(num_obs);
    }
    
    if (y_values.contains(var_name) && dy_values.contains(var_name)) {
        auto y_array = y_values[var_name].get<std::vector<double>>();
        auto dy_array = dy_values[var_name].get<std::vector<double>>();
        for (size_t j = 0; j < num_obs; j++) {
            y_all[j].push_back(y_array[j]);
            dy_all[j].push_back(dy_array[j]);
        }
    } else {
        // Missing observation - fill with NaN
        for (size_t j = 0; j < num_obs; j++) {
            y_all[j].push_back(NaN_VALUE);
            dy_all[j].push_back(NaN_VALUE);
        }
    }
}
```

### 2. **Residual Computation (Model blocks)**

The residual computation in `update_gradient()` methods (e.g., `BloodVessel::update_gradient`, `HybridJunction::update_gradient`) accesses observations like:

```cpp
auto y0 = y[global_var_ids[0]];  // Direct array access
```

**Modification needed:**
```cpp
// Check if observation is valid before using
auto y0_val = y[global_var_ids[0]];
if (std::isnan(y0_val)) {
    // Skip this residual contribution or set to zero
    return;  // Or continue to next variable
}
auto y0 = y0_val;
```

**Better approach:** Modify the residual computation to only compute residuals for equations that involve observed variables. This requires:

1. **Track which variables are observed** (pass `has_observation` vector to `update_gradient`)
2. **Skip residual contributions** for unobserved variables
3. **Adjust Jacobian size** accordingly

### 3. **Optimizer Modifications (LevenbergMarquardtOptimizer.cpp)**

The optimizer currently assumes:
- `y_obs[i]` has size = `num_vars` (all variables)
- `dy_obs[i]` has size = `num_vars` (all variables)

**Option A: Use NaN placeholders**
- Keep current structure
- Check for NaN in residual computation
- Skip residual contributions involving NaN values

**Option B: Sparse observation structure**
- Only include observed variables in `y_obs`/`dy_obs`
- Create mapping: `variable_index -> observation_index`
- Modify residual computation to use mapping

### 4. **Recommended Implementation**

**Step 1: Modify observation reading** to use NaN placeholders (simplest)

**Step 2: Add helper function** to check if variable is observed:
```cpp
bool is_observed(const std::vector<double>& y_vec, size_t var_idx) {
    return !std::isnan(y_vec[var_idx]);
}
```

**Step 3: Modify `update_gradient()` in each block type** to check for NaN:
```cpp
void BloodVessel::update_gradient(...) {
    // Check if all required variables are observed
    bool all_observed = true;
    for (size_t idx : global_var_ids) {
        if (std::isnan(y[idx]) || std::isnan(dy[idx])) {
            all_observed = false;
            break;
        }
    }
    
    if (!all_observed) {
        // Skip residual computation for this block
        // Or compute partial residual using only observed variables
        return;
    }
    
    // ... existing residual computation
}
```

**Step 4: Adjust residual vector size** - The residual vector size is `num_obs * num_eqns`. If we skip some equations, we need to:
- Track which equations contribute to residual
- Adjust Jacobian size accordingly
- Or: Keep full size but set skipped residuals to zero

### 5. **Challenges**

1. **Equation coupling**: Some equations involve multiple variables. If one variable is missing, how do we handle the equation?
   - **Solution**: Skip the entire equation if any required variable is missing
   - **Alternative**: Use estimated/interpolated values for missing variables

2. **Jacobian size**: If we skip equations, the Jacobian size changes
   - **Solution**: Keep full size, set skipped rows to zero
   - **Alternative**: Use sparse structure with only active equations

3. **Convergence**: Fewer observations may affect optimization convergence
   - **Solution**: Adjust tolerances or add regularization

### 6. **Simplest Implementation (Recommended)**

1. **Use NaN placeholders** for missing observations
2. **Skip residual computation** for blocks with missing observations
3. **Keep full residual/Jacobian size** but set skipped contributions to zero
4. **Add warning messages** when observations are missing

This minimizes changes to existing code while enabling optional observations.

