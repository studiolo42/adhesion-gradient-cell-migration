# Adhesion-gradient cell migration model: simulation code

MATLAB code accompanying the manuscript **[TODO: manuscript title]** ([TODO: authors], [TODO: journal, year]).

The script simulates a two-ended (front/back) cell moving along an adhesion gradient and produces the time series shown in the multi-panel figure of the manuscript (position, local flow, and proportion of bound linkers vs. time) for different values of the base parameter β₀.

## Contents

| File | Description |
|------|-------------|
| `adhesion_gradient_model.m` | Complete simulation, data export and plotting script |
| `README.md` | This file |

## Requirements

- MATLAB R2018b or later (uses `yline`); local functions inside the script require R2016b or later.
- Base MATLAB only; no additional toolboxes are used.
- Tested with: [TODO: MATLAB version / OS].

## How to run

1. Open `adhesion_gradient_model.m`.
2. In **Section 1 (User settings)** set `b0` (β₀) and the total simulated time `T`.
3. Run the script. All outputs are written to an `output/` folder in the current working directory.

To reproduce the four columns of the figure, run the script once for each setting below (all other parameters unchanged):

| Panel | β₀ (`b0`) | Approx. simulated time `T` |
|-------|-----------|----------------------------|
| a | 2   | 600 |
| b | 2.5 | 310 |
| c | 6.5 | 72  |
| d | 9   | 21  |

> **[TODO: author to verify]** The `T` values above were read from the time axes of the figure. Please confirm them, and confirm that all other parameters are identical across panels. The script as provided is set to panel d (`b0 = 9`, `T = 20`).

## Model (as implemented)

State variables for the front (f) and back (b) ends:

| Symbol | Code | Meaning |
|--------|------|---------|
| x_f, x_b | `xf_vec`, `xb_vec` | Positions of the front and back |
| x_cm | `x_cm` | Centre-of-mass position, (x_f + x_b)/2 |
| v_f, v_b | `vf_vec`, `vb_vec` | Local flow at the front and back |
| n_f, n_b | `nf_vec`, `nb_vec` | Proportion of bound linkers |
| r_f, r_b | `rf_vec`, `rb_vec` | Adhesion value at each end |

With l = x_f − x_b, Δv = v_f − v_b:

```
r(x)     = A · exp(−K_g · |x − x_peak|)
β_i      = β₀ · (1 + α · r_i / (r_i + r₀))

dx_f/dt  = [ r_f/(r_f + r₀) · v_f − k (l − 1) ] / Γ_f
dx_b/dt  = [ −r_b/(r_b + r₀) · v_b + k (l − 1) ] / Γ_b

dn_f/dt  = r_f (1 − n_f) − n_f · exp{ [ −r_f/(r_f + r₀) · v_f + k (l − 1) ] / (f_s · n_f) }
dn_b/dt  = r_b (1 − n_b) − n_b · exp{ [ −r_b/(r_b + r₀) · v_b + k (l − 1) ] / (f_s · n_b) }

dv_f/dt  = −δ [ v_f − β_f / (1 + (c Δv / d) / (exp(Δv l / d) − 1)) ]  + σ_v · ξ
dv_b/dt  = −δ [ v_b − β_b / (1 + (c Δv / d) / (1 − exp(−Δv l / d))) ] + σ_v · ξ
```

ξ is Gaussian white noise, integrated with the Euler–Maruyama scheme (noise term `sig_v * randn * sqrt(delta*dt)` per step). When v_f = v_b exactly, the 0/0 expression is replaced by its analytical limit, l/(c + l).

The friction factors Γ_f and Γ_b are implemented in the local functions `calculate_Gammaf` and `calculate_Gammab` at the end of the script. For an advancing end, Γ = r/(r + r₀); for a retracting end, Γ = n · k_p · exp(·) (see code for the exact expression).

The simulation stops when either end reaches `x_end`, or when `T` is reached.

## Parameters

| Code name | Value | Description |
|-----------|-------|-------------|
| `b0` | set per run (2, 2.5, 6.5, 9) | Base β (β₀) |
| `c`, `d` | 3.95, 3.95 | Linker-related constants in the flow equations |
| `delta` | 200 | Relaxation rate of the local flows (δ) |
| `k` | 0.8 | Spring constant between front and back |
| `fs` | 5 | Force scale of linker unbinding (f_s) |
| `kp` | 20 | Prefactor of the retraction friction (k_p) |
| `sig_v` | 0.01 | Noise amplitude of the local flows (σ_v) |
| `n_min` | 1e-5 | Lower bound on n (keeps n > 0) |
| `r0` | 0.75 | Half-saturation constant, r/(r + r₀) |
| `alpha` | 0.65 | Strength of β's dependence on r (α) |
| `A` | 27 | Peak value of r |
| `K_g` | 0.1769 | Steepness of the adhesion gradient |
| `cell_track` | 65 | Domain length; `x_peak = cell_track/2`, `x_end = 2·x_peak` |
| `dt` | 1e-4 | Time step |

Initial conditions: x_b = 15, x_f = l_c + 20 with l_c = c / (√(c·β₀/d) − 1), n_f = n_b = 0.25, v_f = v_b = β₀ l / (c + l).

## Outputs

All files are written to `output/`; `<tag>` is `beta0_<value>` (e.g. `beta0_9`, `beta0_2p5`).

| File | Content |
|------|---------|
| `results_<tag>.mat` | Full-resolution data and a `params` struct with every parameter used |
| `results_<tag>.csv` | Down-sampled time series (every `save_stride`-th point). Columns: `t, xf, xb, x_cm, l, vf, vb, nf, nb, rf, rb` |
| `position_vs_time_<tag>.{jpg,svg,fig}` | Position [a.u] vs. Time [a.u] (x_f, x_cm, x_b) |
| `local_flow_vs_time_<tag>.{jpg,svg,fig}` | Local flow [a.u] vs. Time [a.u] (v_f, v_b) |
| `bound_linkers_vs_time_<tag>.{jpg,svg,fig}` | Prop. of bound linkers vs. Time [a.u] (n_f, n_b) |

All quantities are in arbitrary units (a.u.). In the position plots, positions are shown relative to `x_peak` (dashed line at 0), as in the manuscript figure; set `plot_relative_to_peak = false` for absolute positions. Note that `rf`, `rb` in the output are stored with a one-step lag relative to the positions (see the loop in Section 5 of the script).

## Reproducibility

The local-flow equations are stochastic. The random number generator is seeded by `rng_seed` (default 1), so a given seed and parameter set always reproduces the same trajectory. Individual trajectories, especially for large β₀, depend on the seed. **[TODO: author to state the seed used for the manuscript figures, or note that figures show a representative run.]**

## Runtime and memory

Runtime and memory scale with `T/dt`. `T = 20` is 2×10⁵ steps and is quick. `T = 600` is 6×10⁶ steps; expect a long run and roughly 0.5 GB of memory for the stored vectors.

## License and citation

[TODO: license, e.g. MIT / CC BY 4.0]

If you use this code, please cite: [TODO: manuscript citation / DOI].

## Contact

[TODO: name, email]
