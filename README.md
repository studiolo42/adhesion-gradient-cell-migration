# Four migration modes: simulation code

MATLAB code for Fig. 3 of the manuscript. It simulates a cell with a front and a back end moving along an adhesion gradient.

## Files
- `adhesion_gradient_model.m`: the simulation and plots
- `README.md`: this file

## How to run
1. Open `adhesion_gradient_model.m` in MATLAB (R2018b or newer, no toolboxes needed).
2. Set `b0` and `T` at the top (see table below).
3. Press Run.

## Settings for the four panels

| Panel | `b0` | `T` |
|-------|------|-----|
| a | 2   | 600 |
| b | 2.5 | 310 |
| c | 6.5 | 72  |
| d | 9   | 21  |

The file is set to panel d. All other parameters stay the same. Large `T` (e.g. 600) takes a while to run.

## What you get
Three figures, all in arbitrary units (a.u):
1. **Position vs. Time**: front (x_f, red), centre (x_cm, green) and back (x_b, blue). The dashed line is the adhesion peak.
2. **Local flow vs. Time**: front (v_f, red) and back (v_b, blue).
3. **Proportion of bound linkers vs. Time**: front (n_f, red) and back (n_b, blue).

## Notes
- The local flow has random noise, so every run looks slightly different.
- The run stops early if the cell reaches the end of the track.

## Contact
For any further queries, please write to: papagantijeevankumar@gmail.com
