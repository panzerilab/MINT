import matlab.engine
import numpy as np

# Start the MATLAB engine
eng = matlab.engine.start_matlab('-nojvm')

# Add the MINT toolbox path to MATLAB (check if this path is correct)
eng.eval('addpath(genpath(fullfile(\'/home\', \'matias\', \'Documents\', \'Projects\', \'MINT\')));', nargout=0)
eng.eval('disp(matlabroot)', nargout=0)

# Check if PID function exists
result = eng.eval('exist("PID", "file")')
if result != 2:
    print("Error: PID function is not available. Check your MINT toolbox setup.")
else:
    print("PID function is available.")

# Initialize inputs for PID function
X = np.random.choice([0, 1], size=(5, 100))  # 5 neurons, 100 trials
Y = np.random.choice([0, 1], size=(5, 100))  # 5 neurons, 100 trials
Z = np.random.choice([0, 1], size=(5, 100))  # 5 neurons, 100 trials

# Convert the data into MATLAB types (matlab.double)
X_matlab = matlab.double(X.tolist())  # Convert NumPy array to MATLAB double
Y_matlab = matlab.double(Y.tolist())  # Convert NumPy array to MATLAB double
Z_matlab = matlab.double(Z.tolist())  # Convert NumPy array to MATLAB double

# Combine X and Y into a list (not numpy array) as required by MATLAB
XYZ = [X_matlab, Y_matlab, Z_matlab]  # List of MATLAB double arrays

# Define the required outputs you want to compute
req_outputs = ['PID_atoms']  # Example outputs

# Convert Python lists to MATLAB cell arrays using eng.cell()
bin_method = eng.cell([['none']])  # MATLAB cell array of 1x1 cell containing 'none'
n_bins = matlab.double([3])  # Convert to MATLAB double array

# Define optional arguments as a dictionary (convert to matlab types where needed)
opts = {
    'bias': 'shuffSub',  # No correction
    'shuff': matlab.double(5),
    'pid_constrained': False,
    'bin_method': bin_method,  # MATLAB cell array
    'n_bins': [n_bins],  # MATLAB double array
    'computeNulldist': True,  # Boolean for null distribution
    'n_samples': matlab.double(5),  # Number of samples for null distribution
    'suppressWarnings': False,  # Boolean to show/hide warnings
    'NaN_handling': 'removeTrial'  # Handle NaN values by removing trials
}

# If PID is available, call the PID function
if result == 2:
    PID_values, PID_naive, PID_nullDist= eng.PID(XYZ, req_outputs, opts, nargout=3) #, PID_naive, PID_nullDist
    # Process the results returned by the PID function
    print("PID Values: ", np.array(PID_values))
    print("Naive PID Estimates: ", np.array(PID_naive))
    print("Null Distribution: ", np.array(PID_nullDist))
else:
    print("PID function not found. Please check your MINT toolbox setup.")

# Close the MATLAB engine
eng.quit()
