# MINT Multivariate Information in Neuroscience Toolbox
## Software requirements
#### Linux:
- MATLAB supported C/C++ compiler (tested with gcc and g++ version > 6.3.1)
- MATLAB (tested with versions > r2018b but earlier version should be compatible)
  - Statistics and Machine Learning toolbox
  - Optimization toolbox
  - Parallel Computing Toolbox
  - Signal Processing Toolbox

#### MacOS:
- MATLAB supported C/C++compiler
- MATLAB (tested with versions > r2018b but earlier version should be compatible)
  - Statistics and Machine Learning toolbox
  - Optimization toolbox
  - Parallel Computing Toolbox
  - Signal Processing Toolbox

#### Windows:
- MATLAB supported C/C++ compiler (tested with MinGW-w64, see installation instructions [here](https://it.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler))
- Microsoft Visual C++
- Microsoft Visual C++ build tools
- MATLAB (tested with versions > r2018b but earlier version should be compatible)
  - Statistics and Machine Learning toolbox
  - Optimization toolbox
  - Parallel Computing Toolbox
  - Signal Processing Toolbox

## Build instructions
1. Before starting the installation make sure that all software requirements described above are satisfied.
2. After this, open MATLAB, navigate to the main MINT directory, and run the script `BuildMINT.m`.
    - The script compiles the code that needs compiling, adds the required source files to MATLAB search path by adding a `StartupMINT.m` file in the MATLAB user home folder and by editing the default user `Startup.m` file and performs the required install tests. If all software requirements are satisfied the process should seamlessly finish without errors.
    
### Build with Docker
1. If you have Docker installed, you can build the Dockerfile in this package by running 
``` 
docker build --build-arg ADDITIONAL_PRODUCTS="Statistics_and_Machine_Learning_Toolbox Optimization_Toolbox Parallel_Computing_Toolbox Signal_processing_Toolbox" -t matlab4mint:R2024b .
```
2. After the Docker image has been built with all the necessary toolboxes, run 
```
docker run --init -v /path/on/host:/home/matlab/shared -e MLM_LICENSE_FILE=27000@MyServerName matlab4mint:R2024b 
``` 
Where /path/on/host can be replaced with the directory in your computer where you can work persistently after the Docker image is used. If you do not know your license hostname, you can skip -e MLM_LICENSE_FILE=27000@MyServerName and you can login with your username and password when you are inside the container.
3. Once inside the container, write
```
cd /home/matlab/MINT
run('BuildMINT.m')
``` 

## License
MINT is licensed under the GNU General Public License (GPL), version 3 or later.
This means that you are free to redistribute and modify the software under the terms of the GNU GPL, either version 3 of the License or (at your option) any later version, as published by the Free Software Foundation.
MINT is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).

### DOI 
[![DOI](https://zenodo.org/badge/872930487.svg)](https://doi.org/10.5281/zenodo.13998526)

## How to Use MINT
To get started with MINT, refer to the tutorials provided in the `How_to_use_MINT` folder. This folder contains a range of step-by-step guides that demonstrate how to use the main functions of the MINT toolbox. 





