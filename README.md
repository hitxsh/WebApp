Make sure you have MATLAB installed in your local device.

## Create Environment

python -m venv myworld


myworld\Scripts\activate

## Install Django

python3.10 -m pip install Django

## Install required packages

python3.10 -m pip uninstall -y numpy scipy open3d trimesh pyglet

python3.10 -m pip install Pillow psycopg2 numpy==1.25.0 scipy==1.10.0 trimesh open3d==0.18.0 pyglet==1.5.27


## Go to Imagedisplay folder:

cd Imagedisplay

## Run the server:

python3.10 manage.py runserver
