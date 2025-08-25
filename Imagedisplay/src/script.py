import matlab.engine
import os
import shutil

# Setup paths
code_folder = '/home/vamshi/Desktop/NGU_WebApp/MatlabCode'
input_path = '/home/vamshi/Desktop/NGU_WebApp/media/input.tif'
output_path = '/home/vamshi/Desktop/NGU_WebApp/media/output.tif'

# Start MATLAB and add the code path
eng = matlab.engine.start_matlab()
eng.addpath(code_folder, nargout=0)

# Run your wrapped script
eng.run_main(input_path, output_path, nargout=0)

# Quit MATLAB
eng.quit()

# Confirm output exists
if os.path.exists(output_path):
    print("✅ MATLAB processing complete. Output saved at:", output_path)
else:
    print("❌ Output not found. Check Codemsir script.")
