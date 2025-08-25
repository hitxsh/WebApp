# import matlab.engine
# import sys

# # Start MATLAB engine
# eng = matlab.engine.start_matlab()

# # Optional: change directory to where your .m file is located
# eng.cd(r'C:\Users\Hitesh\Downloads\testv1\project')

# # Run the MATLAB file
# eng.workspace['filename'] = sys.argv[1]
# print(f"Running MATLAB script with filename: {sys.argv[1]}")
# eng.Codemsir(nargout=0)  # Use nargout=0 if script doesn't return anything

# # Stop MATLAB
# eng.quit()

import os, sys, shutil
src_path = sys.argv[1]
print(f"Running MATLAB script with filename: {src_path}")

# Change into the project directory
project_dir = r"./src"
filename = os.path.basename(src_path)
shutil.copy(src_path, project_dir)
os.chdir(project_dir)
os.system(f"matlab -batch \"Codemsir('{filename}');\"")