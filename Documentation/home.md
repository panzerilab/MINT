# &#128596; MINT Multivariate Information in Neuroscience Toolbox
## Software requirements
#### Linux:
- MATLAB supported C/C++ compiler (tested with gcc and g++ version > 6.3.1)
- MATLAB (tested with versions > r2019b but earlier version should be compatible)
  - Statistics and Machine Learning toolbox
  - Optimization toolbox
  - Parallel Computing Toolbox
  - Communications Toolbox
  - DSP System Toolbox

#### MacOS:
- MATLAB supported C/C++compiler
- MATLAB (tested with versions > r2019b but earlier version should be compatible)
  - Statistics and Machine Learning toolbox
  - Optimization toolbox
  - Parallel Computing Toolbox
  - Communications Toolbox
  - DSP System Toolbox

#### Windows:
- MATLAB supported C/C++ compiler (tested with MinGW-w64, see installation instructions [here](https://it.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler))
- Microsoft Visual C++
- Microsoft Visual C++ build tools
- MATLAB (tested with versions > r2019b but earlier version should be compatible)
  - Statistics and Machine Learning toolbox
  - Optimization toolbox
  - Parallel Computing Toolbox
  - Communications Toolbox
  - DSP System Toolbox

## Build and test instructions
1. Before starting the installation make sure that all software requirements described above are satisfied.
2. After this open MATLAB and run the script `BuildAndTest.m`

The script compiles the code that needs compiling, adds the required source files to MATLAB search path by adding a `StartupNIT.m` file in the MATLAB user home folder and by editing the default user `Startup.m` file and performs the required install tests. If all software requirements are satisfied the process should seamlessly finish without errors.

## License
