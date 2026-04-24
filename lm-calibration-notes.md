# Levenberg–Marquardt calibration: when residual stays high

## What the current LM is (and is not) optimizing for

In `LevenbergMarquardtOptimizer::run`, iterations stop when **both**:

- \(\|J^\top r\|\) (logged as **norm grad**), and
- \(\|\Delta\alpha\|\) (logged as **norm inc**)

are below tolerances (`src/optimize/LevenbergMarquardtOptimizer.cpp`).

That is a **stationarity / “no progress”** criterion, **not** “\(\|r\|\) is small.”

You can therefore see **small `norm grad` / `norm inc`** but **large `norm residual`**: flat directions, poor scaling, rank deficiency, or **no** parameter vector \(\alpha\) that fits the data with the chosen circuit model.

**Suggestion:** treat **fit quality** explicitly—for example monitor **`norm residual`** (already printed) or cost \(\tfrac{1}{2}\|r\|^2\)—separately from the current stopping rule, and optionally add a **residual-based** stop or warning.

---

## Implementation-specific notes (this codebase)

### 1. Full step only (no line search)

After `delta = mat.llt().solve(vec)` the code does `alpha -= delta` with **no** backtracking or accept/reject.

If the linearization is poor, steps can **overshoot** or **oscillate** while \(\|r\|\) remains large.

**Mitigations:** line search on \(\|r(\alpha - t\Delta)\|\), or trust-region LM that only accepts steps that reduce the model.

### 2. \(\lambda\) update heuristic

\(\lambda\) is updated as `lambda *= vec.norm() / vec_old.norm()` with `vec = J^T r`.

That is a **simple ratio heuristic**; it can change \(\lambda\) quickly and is **not** tied to whether \(\|r\|\) actually decreased.

**Mitigation:** classical LM schemes adjust \(\lambda\) from **step success** (increase if residual goes up, decrease if down), for example Nielsen-style rules.

### 3. Normal equations and `LLT`

The solve uses `LLT` on \(J^\top J + \lambda \operatorname{diag}(J^\top J)\).

If \(J\) is **rank-deficient** or **badly scaled**, the normal matrix is **ill-conditioned**; the factorization can **fail** or yield poor \(\Delta\alpha\).

**Mitigations:** LDLT with pivoting, SVD-based pseudo-inverse, or **column scaling** of \(J\).

### 4. Fixed parameters (`fixed_param_ids`)

The code solves the full system then **zeros selected components of `delta`**.

That is **not** equivalent to LM in the **free-parameter subspace** and can **distort** the step for remaining parameters.

**Mitigation:** build \(J\) with **only free columns**, or impose constraints in a KKT-style system.

### 5. `NORMAL_JUNCTION` in calibration

`Junction::update_gradient` contributes **two** residual rows and **no** Jacobian columns (`src/model/Junction.cpp`).

That can change the **effective** least-squares problem versus the time-domain DAE assembly.

If the network is junction-heavy, **calibration versus forward-model mismatch** can contribute to a floor on \(\|r\|\).

---

## Problem / data side (often the real limiter)

- **Inconsistent or noisy** \((y, \dot y)\) versus a **pure** 0D element model.
- **Units / time base** errors, or **partial / NaN** observations (some blocks skip residuals when data are NaN).
- **Too little freedom** in \(\alpha\) (for example aggressive **`freeze_connector_segments`**, fixed \(C\) via **`BloodVesselFC`**, and so on).
- **Non-identifiability**: many \(\alpha\) give similar \(y\) at your sample times, leading to small \(\|J^\top r\|\) but **nonzero** \(\|r\|\).

---

## Cheap knobs to try first

| Idea | Direction |
|------|-----------|
| `initial_damping_factor` (\(\lambda_0\)) | Increase if unstable or huge steps; decrease if barely moving (carefully). |
| `maximum_iterations` | Raise if stopping on tolerance while \(\|r\|\) is still large. |
| `tolerance_gradient`, `tolerance_increment` | Tighten if you stop “successfully” but fit is poor; loosen if never finishing. |
| **Initial \(\alpha\)** | Often more impactful than LM tweaks alone. |

---

## Stronger changes (more engineering)

- **Row-wise weighting** of residuals so pressure versus flow contribute comparably to \(\|r\|\).
- **Nielsen-style \(\lambda\)** with **gain ratio** \(\rho\) and step acceptance.
- **Box constraints** (for example \(C > 0\)): consider **L-BFGS-B** or projected steps if bounds dominate the failure mode.

---

## Diagnostic split (high level)

- If **`norm grad` / `norm inc` are small** but **`norm residual` is large** → think **scaling, identifiability, model/data mismatch, fixed-parameter projection**, or **wrong stopping criterion**.
- If **`norm grad` stays large** → **ill-conditioning**, **bad steps**, or **inconsistent Jacobian versus residual** definition.

---

## File references

- `src/optimize/LevenbergMarquardtOptimizer.cpp` — iteration loop, stopping rule, \(\lambda\) update, `llt` solve, fixed `delta` zeroing.
- `src/optimize/LevenbergMarquardtOptimizer.h` — mathematical formulation in comments.
- `src/model/Junction.cpp` — `NORMAL_JUNCTION` LM residual contribution (`update_gradient`).
