# Computational Structure of Egocentric Boundary Cell Responses in Retrosplenial Cortex

This repository contains the MATLAB implementation and demo data for the computational modeling of **Egocentric Boundary Cells (EBCs)** in the mouse **Retrosplenial Cortex (RSC)**. The demo is designed to be lightweight and to reproduce **Figure 1** from the manuscript by running a single script.


## Quick start (reproduce Figure 1)

1. Open MATLAB.
2. Set your MATLAB current folder to the repository root.
3. Run:

```matlab
   run('Demo.m')
```
> **Note:** The demo is intended to run on a typical desktop/laptop. Full-scale fitting across many cells/sessions may require cluster computing.

---

## Data Format
The demo includes sample data files with the following structures:

### 1. `spike_timestamps.txt`
* **Content**: A single-column list of timestamps representing each action potential (spike).
* **Units**: Microseconds ($\mu s$).

### 2. `tracking_data.txt`
* **Content**: Video tracking data with one row per video frame.
* **Columns**: `[timestamp, x_position, y_position, head_direction]`
    * **Timestamp**: Microseconds ($\mu s$).
    * **X / Y Position**: Animal's coordinates in centimeters (cm).
    * **Head Direction**: Degrees (increasing counterclockwise from the positive x-axis).

---

## Requirements
* **MATLAB** (Recommended: R2021a or later)
* **Optimization Toolbox**: Required for fmincon.
* **Global Optimization Toolbox**: Required for GA.
* **Parallel Computing Toolbox**: Required for parallelizing model fitting across multiple repetitions (`nworkers`).

---

## Configuration
You can adjust the following parameters within `Demo.m` to fit your computational resources or research needs:

| Parameter | Description | Default |
| :--- | :--- | :--- |
| `rf` | List of receptive field models to fit (e.g., `dog`, `bounded_gaussian`). | `{'angle_distance_gaussian', ...}` |
| `bps` | Number of boundary points included (1 to 120). | `1` |
| `r` | Number of fitting repetitions per model for stability. | `5` |
| `nworkers` | Number of CPU cores for parallel processing. | `4` |
| `fps` | Frame rate of tracking data for firing rate visualization. | `30` |

---

## Outputs
The script automatically organizes results into the following directory structure:

* **`/Demo results/`**: Contains `.mat` files for each repetition and the final `results_best.mat` (selected via log-likelihood optimization).
* **`/Demo figures/`**: Contains generated visualizations.
    * **`/{bps}pts/`**: Subfolders sorted by the number of boundary points used.
    * **`..._rf.png`**: The final receptive field visualizations.

