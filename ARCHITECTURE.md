# HiFAST Code Architecture

The source code structure of the HiFAST project is designed with clear modularity, separating the data processing stages, core algorithms, utilities, and interactive interfaces.

The structure can be summarized in the following main components:

### 1. Top-Level Pipeline Scripts

*   **Location:** `hifast/*.py` (e.g., `flux.py`, `bld.py`, `rfi.py`, `sw.py`, `cube.py`)
*   **Role:** These scripts serve as the **main entry points** and drivers for the HiFAST data processing pipeline. Each file corresponds to a distinct, complete processing stage (e.g., flux calibration, baseline fitting, cube synthesis). Users invoke these via the command line (e.g., `python -m hifast.flux ...`), which makes the pipeline staged and user-friendly.

### 2. Core Algorithm Library (`core/`)

*   **Location:** `hifast/core/`
*   **Role:** This is the **"engine room"** of the project, containing all the core, reusable scientific computing and astronomical algorithms. For example, `baseline.py` implements baseline fitting algorithms, `grid.py` handles data gridding, and `cal.py` contains universal calibration logic. The top-level scripts call functions from this library to perform the actual computations. This design cleanly separates "flow control" from "core computation."

### 3. Specialized Modules

*   **Location:** `hifast/cbr/` and `hifast/ripple/`
*   **Role:** These are specialized modules created to solve specific, complex problems.
    *   `cbr/`: Focuses on the logic for **Calibrator**-based flux scaling, providing critical support for `flux.py`.
    *   `ripple/`: Dedicated to handling **Standing Waves** and **RFI**, a key part of the data cleaning process.
*   These modules also demonstrate the project's capacity to manage and integrate solutions for distinct sub-problems.

### 4. Utility Library (`utils/`)

*   **Location:** `hifast/utils/`
*   **Role:** This directory houses general-purpose helper tools that are not tied to a specific scientific algorithm. For instance, `io.py` manages file I/O, `conf_arg.py` handles command-line argument parsing, and `output.py` manages logging. This promotes code reuse and improves development efficiency.

### 5. Interactive Modules (`interaction/`)

*   **Location:** `hifast/interaction/`
*   **Role:** A standout feature of the project, this directory is dedicated to enabling **interactive data processing** within a Jupyter Notebook environment. Files like `widgets.py` and `bld_i.py` indicate the creation of graphical user interfaces for debugging complex parameters, which significantly enhances usability.

### 6. Automation and Batch-Processing Scripts (`scripts/`)

*   **Location:** `scripts/`
*   **Role:** This directory contains high-level wrapper scripts for automating and orchestrating the pipeline. The primary script, `hifast.sh`, acts as a powerful batch-processing tool. It allows users to run a sequence of `hifast` commands (defined in a parameter file) on multiple data files in parallel, managing the workflow and chaining the output of one step to the input of the next. This provides a crucial layer of automation for large-scale survey processing.