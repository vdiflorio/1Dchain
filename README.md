# 1Dchain

C++ code for the numerical simulation of one-dimensional oscillator chains and
the computation of transport properties using the **Transient Time Correlation Function (TTCF)** formalism.

In addition, the code computes the **direct time average of the same observables**
used in the TTCF, allowing for a systematic comparison between TTCF predictions and time-averaged quantities.

Post-processing and parameter handling are performed using Python scripts.

---

## Directory structure

```directory
1Dchain/
├── src/                # C++ source files
│   ├── main.cpp        # Main executable
│   ├── ode_func.cpp    # Equations of motion
│   ├── ode_solvers.cpp # Time integration schemes
│   ├── param.cpp       # Simulation parameter handling
│   ├── utils.cpp       # Utility functions
│   └── Makefile        # Build system
├── include/            # C++ header files
│   ├── ode_func.h
│   ├── ode_solvers.h
│   ├── param.h
│   └── utils.h
├── examples/           # Example input files
│   └── parametri.json  # Sample simulation parameters
├── external/           # External dependencies
│   └── json.hpp        # nlohmann/json single-header library
├── python_script/      # Python post-processing scripts
│   ├── change_parameter.py
│   └── read_TTCF_results.ipynb
└── README.md
```

---

## Compilation

The code is written in **C++17** and is compiled using a simple `Makefile`.

From the `src/` directory, run:

```bash
make clean all
```

This produces the executable in the same directory.

---

## Running a simulation

From the src/ directory:

```bash
./fput ../examples/parametri.json
```

The input file specifies:
 - the physical parameters of the 1D oscillator system,
 - numerical integration settings,
 - the observables used for TTCF and time-average calculations.

---

## Authors
 - Vincenzo Di Florio
 - Davide Carbone

---

## License

This project is licensed under the GNU General Public License v3.0 or later.

---
