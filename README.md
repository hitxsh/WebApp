# Design_credit
Create Environment:
python -m venv myworld


myworld\Scripts\activate.bat

Install Django:

python -m pip install Django

Install required packages

python -m pip uninstall -y numpy scipy open3d trimesh pyglet
python -m pip install Pillow psycopg2 numpy==1.25.0 scipy==1.10.0 trimesh[all] open3d==0.18.0 pyglet==1.5.27 numpy


Go to Imagedisplay folder & run:

python manage.py runserver


Check this:
https://www.w3schools.com/django/django_getstarted.php
