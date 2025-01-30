import matlab.engine
import numpy as np

# Start the MATLAB engine
eng = matlab.engine.start_matlab('-nojvm')

# Add the MINT toolbox path to MATLAB (check if this path is correct)
eng.eval('addpath(genpath(fullfile(\'/home\', \'username\', \'Documents\', \'Projects\', \'MINT\')));', nargout=0)
eng.eval('disp(matlabroot)', nargout=0)

# Check if MI function exists
result = eng.eval('exist("MI", "file")')
if result != 2:
    print("Error: MI function is not available. Check your MINT toolbox setup.")
else:
    print("MI function is available.")

# Initialize inputs for MI function
X = np.random.choice([0, 1], size=(10, 100))  # 10 neurons, 100 trials
Y = np.random.choice([0, 1], size=(10, 100))  # 10 neurons, 100 trials

# Convert the data into MATLAB types (matlab.double)
X_matlab = matlab.double(X.tolist())  # Convert NumPy array to MATLAB double
Y_matlab = matlab.double(Y.tolist())  # Convert NumPy array to MATLAB double

# Combine X and Y into a list (not numpy array) as required by MATLAB
XY = [X_matlab, Y_matlab]  # List of MATLAB double arrays

# Define the required outputs you want to compute
req_outputs = ['I(A;B)', 'Ilin(A;B)', 'coI(A;B)']  # Example outputs

# Convert Python lists to MATLAB cell arrays using eng.cell()
bin_method = eng.cell([['none']])  # MATLAB cell array of 1x1 cell containing 'none'
n_bins = matlab.double([3])  # Convert to MATLAB double array

# Define optional arguments as a dictionary (convert to matlab types where needed)
opts = {
    'bias': 'pt',  # No correction
    'bin_method': bin_method,  # MATLAB cell array
    'n_bins': [n_bins],  # MATLAB double array
    'computeNulldist': True,  # Boolean for null distribution
    'n_samples': 10,  # Number of samples for null distribution
    'suppressWarnings': False,  # Boolean to show/hide warnings
    'NaN_handling': 'removeTrial'  # Handle NaN values by removing trials
}

# If MI is available, call the MI function
if result == 2:
    MI_values, MI_naive, MI_nullDist = eng.MI(XY, req_outputs, opts, nargout=3) #
    # Process the results returned by the MI function
    print("MI Values: ", np.array(MI_values))
    print("Naive MI Estimates: ", np.array(MI_naive))
    print("Null Distribution: ", np.array(MI_nullDist))
else:
    print("MI function not found. Please check your MINT toolbox setup.")

# Close the MATLAB engine
eng.quit()
