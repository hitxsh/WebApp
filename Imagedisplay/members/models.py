from django.db import models
from django.contrib.auth.models import User

class Profile(models.Model):
    user = models.OneToOneField(User, on_delete=models.CASCADE)
    firstname = models.CharField(max_length=255)
    lastname = models.CharField(max_length=255)

    def __str__(self):
        return self.user.username

class UploadedImage(models.Model):
    image = models.ImageField(upload_to='uploads/')
    uploaded_at = models.DateTimeField(auto_now_add=True)


class UploadedSTL(models.Model):
    stl_file = models.FileField(upload_to='stl_uploads\\')
    uploaded_at = models.DateTimeField(auto_now_add=True)

# models.py
# from django.db import models

# class STLFile(models.Model):
#     stl_file = models.FileField(upload_to='stl_files/')
#     uploaded_at = models.DateTimeField(auto_now_add=True)
    
#     def __str__(self):
#         return self.stl_file.name

# In your models.py

import os
from django.contrib.auth.models import User
from django.db import models
from django.contrib.postgres.fields import JSONField  # Or models.JSONField in Django 3.1+

class STLFile(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    run_folder = models.CharField(max_length=512)
    input_file = models.CharField(max_length=512)
    output_file = models.CharField(max_length=512)
    metrics = models.JSONField()  # Store metrics as a dict
    timestamp = models.DateTimeField(auto_now_add=True)

