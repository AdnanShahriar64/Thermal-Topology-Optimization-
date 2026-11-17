# Ground IR + albedo: pre-solver-only MATLAB patch

## Scope

The public repository has one MATLAB script, `Runner_4.m`, containing the
`Thermal_Model_Runner_Complete` local function. The patch inserts a block just
after `K_Conv` is assembled, before the existing surface-normal and solver code.
The original `Surface_Integration_3D_Flux_Thermal` and all lines from
`%% Solver` onward remain **unchanged**.

This delivery intentionally **precomputes but does not yet apply** the two new
loadings to the temperature integration. To activate them, a later change to
`f = fsl + frd + fgr` inside the solver is necessary. That change has **not**
been made because the request expressly excludes solver edits.

## Installation

1. Copy these three helpers next to the original `Runner_4.m`:
   `Ground_Radiation_Precompute.m`, `Ground_Radiation_IR_Load.m`,
   `Ground_Radiation_Albedo_Load.m`.
2. Run:
   `python apply_ground_radiation_patch.py /path/to/Runner_4.m`
3. The script creates `Runner_4_ground_radiation.m` beside the original.
   Run it with the existing Element_Library and Abaqus .inp dependencies.
4. Alternatively, insert the block in
   `Thermal_Model_Runner_Complete_presolver_insert.txt` immediately following
   `[K_Conv] = Assembly(Ke_Convection,1,Mesh_Information,Element_Information);`.

If your function is stored as `Thermal_Model_Runner_Complete.m` instead,
pass that file to the patch script. Its original name remains unchanged.

## New quantities

- `Ground_Radiation.F_Area` [m^2/node] = integral of N_i * F_ground over **outer group 1 only**.
- `f_ground_IR` [W/node] = incoming longwave thermal load, evaluated once using
  the placeholder constant `Ground_Radiation.T_ground_IR = T_ground`.
- `f_albedo` = function handle. At any solar time evaluate
  `f_albedo_now = f_albedo(SF,Sol_Vec)` to obtain a nodal vector [W].
  `Sol_Vec` in the existing code points in the direction sunlight **travels**,
  so solar zenith is computed using `-Sol_Vec`.
- `Ground_Radiation.Info` provides area and view-factor diagnostics.

Current parameters: ground IR emissivity = 0.96, diffuse albedo = 0.12;
wall IR emissivity and solar absorptivity are read from existing material data.

## Surface integration method

Uses existing per-Gauss-point integration outputs:

`T_surf_gp(:,gp,face)` = surface Gauss weight x shape functions x dA;
`dTeS_gp(:,gp,face)` = oriented vector area. `Ke_surf` is not modified.
View factor is computed per Gauss point, which is more accurate than using
one nodal value with the already-integrated `Ke_surf` matrix on a curved face.

`F_ground = (1 - n_z)/2` applies to a convex dome in vacuum above a flat,
unobstructed, infinite horizontal ground plane. The outward normal is selected
using the radial vector from [0,0,0]. Modify the dome center if needed.
If terrain occlusion, self-shadowing, diffuse interreflection, a sloped ground
surface, or a general nonconvex shape matters, this method is insufficient.

## Mathematical model

f_IR = eps_wall * eps_ground * sigma * T_ground_IR^4 * F_Area

f_albedo = alpha_wall * rho_ground * SF * max(0,-Sol_Vec dot up) * F_Area

Both are incident **positive** heating loads. IR wall emission already exists
in the model and is separate from these incoming components. `T_ground_IR`
should eventually be a time-dependent ground temperature, independent from
the fixed-temperature conductive base boundary condition.

## Validation

Run `test_Ground_Radiation` in MATLAB. It uses synthetic surface integration
arrays and checks orientation, surface selection, energy balance, and day/night
albedo. The repository does not currently include its `Element_Library` or
Abaqus mesh input, so it has **not** been possible to run an end-to-end model.
MATLAB/Octave is not available in the artifact environment either.
