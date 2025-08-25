from django.shortcuts import render, redirect
from django.contrib.auth import authenticate, login, logout
from .forms import SignUpForm
from .forms import ImageUploadForm
from .models import UploadedImage
from PIL import Image
import os
from django.conf import settings
from django.contrib.auth.decorators import login_required
from .image_processor import get_image_dimensions
import shutil

@login_required
def home(request):
    if request.user.is_authenticated:
        return render(request, 'home.html', {'username': request.user.username})
    return redirect('login')

def login_view(request):
    if request.method == 'POST':
        username = request.POST['username']
        password = request.POST['password']
        user = authenticate(request, username=username, password=password)
        if user is not None:
            login(request, user)
            return redirect('home')
        else:
            return render(request, 'login.html', {'error': 'Invalid Credentials'})
    return render(request, 'login.html')

def logout_view(request):
    logout(request)
    return redirect('login')

def signup_view(request):
    if request.method == 'POST':
        form = SignUpForm(request.POST)
        if form.is_valid():
            user = form.save()
            return redirect('login')
    else:
        form = SignUpForm()
    return render(request, 'signup.html', {'form': form})



import os
import shutil
from django.conf import settings
from django.shortcuts import render, redirect
from .forms import ImageUploadForm
from .models import UploadedImage
from .image_processor import get_image_dimensions

@login_required
def upload_image(request):
    if request.method == 'POST':
        print("POST request received")
        username = request.user.username
        os.makedirs(os.path.join(settings.MEDIA_ROOT, username), exist_ok=True)
        form = ImageUploadForm(request.POST, request.FILES)
        if form.is_valid():
            # assign the logged-in user before saving
            uploaded_image = form.save()
            print(settings.MEDIA_ROOT)
            original_path = os.path.join(settings.MEDIA_ROOT, str(uploaded_image.image))
            filename = os.path.basename(original_path)
            print(f"Original path: {original_path}")
            
            copy_destination_dir = os.path.join(settings.MEDIA_ROOT, '..\\src')
            copy_destination_dir2 = os.path.join(settings.MEDIA_ROOT, 'processed', username, filename.split('.')[0])
            os.makedirs(copy_destination_dir2, exist_ok=True)
            os.makedirs(copy_destination_dir, exist_ok=True)
            shutil.copy(original_path, copy_destination_dir2)

            exit_status = os.system(f"python3.10 ./src/python_wrapper.py {original_path}")

            width, height = get_image_dimensions(str(uploaded_image.image))
            new_filename = f"output2_{filename}.jpg"
            results_filename = f"results_{filename}.txt"
            last_image_path = f"{copy_destination_dir}\\{new_filename}"
            results_path = f"{copy_destination_dir}\\{results_filename}"

            destination_dir = os.path.join(settings.MEDIA_ROOT, 'processed')
            os.makedirs(destination_dir, exist_ok=True)
            destination_path = os.path.join(destination_dir, new_filename)
            print(last_image_path)
            shutil.copy(last_image_path, destination_path)
            shutil.copy(last_image_path, os.path.join(copy_destination_dir2, new_filename))
            shutil.copy(results_path, os.path.join(copy_destination_dir2, results_filename))
            processed_image_url = f"{settings.MEDIA_URL}processed/{request.user.username}/{filename.split('.')[0]}/{new_filename}"
            results_url = f"{settings.MEDIA_URL}processed/{request.user.username}/{filename.split('.')[0]}/{results_filename}"
            print(results_url)
            with open('.'+results_url, 'r') as f:
                lines = f.read().splitlines()

            print(processed_image_url)

            return render(request, 'upload.html', {
                'form': form,
                'processed_image_url': processed_image_url,
                'roughness': lines[0],
                'roundness': lines[1],
                'sphericity': lines[2],
                'nrq': lines[3]
            })
    else:
        form = ImageUploadForm()
    
    return render(request, 'upload.html', {'form': form})



import os
from django.conf import settings
from django.shortcuts import render
from django.contrib.auth.decorators import login_required

import os
from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.shortcuts import render


@login_required
def recent_results(request):
    user_dir = os.path.join(settings.MEDIA_ROOT, 'processed', request.user.username)
    runs = []

    if os.path.exists(user_dir):
        for folder_name in os.listdir(user_dir):
            folder_path = os.path.join(user_dir, folder_name)
            if os.path.isdir(folder_path):
                result_file = ''
                image_file = ''
                for f in os.listdir(folder_path):
                    if f.startswith('results_') and f.endswith('.txt'):
                        result_file = os.path.join(folder_path, f)
                    if f.startswith('output2_') and (f.endswith('.jpg') or f.endswith('.png')):
                        image_file = os.path.join(folder_path, f)
                    if f.endswith('.tif') or f.endswith('.tiff'):
                        original_image_file = os.path.join(folder_path, f)

                # Skip if necessary files are missing
                if not result_file or not image_file:
                    continue

                # Read results from file
                try:
                    with open(result_file, 'r') as rf:
                        lines = rf.read().splitlines()
                        if len(lines) < 4:
                            continue
                        result_data = {
                            'run_name': folder_name,
                            'original_image_url': os.path.join(settings.MEDIA_URL, 'processed', request.user.username, folder_name, os.path.basename(original_image_file)),
                            'image_url': os.path.join(settings.MEDIA_URL, 'processed', request.user.username, folder_name, os.path.basename(image_file)),
                            'results': {
                                'roughness': lines[0],
                                'roundness': lines[1],
                                'sphericity': lines[2],
                                'nrq': lines[3],
                            }
                        }
                        runs.append(result_data)
                except Exception as e:
                    print(f"Failed to read result for {folder_name}: {e}")
                    continue
    runs.reverse()
    return render(request, 'recent.html', {'runs': runs})


import os
import uuid
import json
import numpy as np
from datetime import datetime
from django.conf import settings
from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from .forms import STLUploadForm, STLProcessingForm
from .models import STLFile
from stl_src.stl_processor import process_stl


def save_upload_log(user, original_filename, save_log_dir):
    log_dir = os.path.join(settings.MEDIA_ROOT, 'stl_processed', f'{user.username}')
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, 'upload_log.txt')

    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(log_file, 'a') as f:
        f.write(f"{original_filename},{save_log_dir},{timestamp}\n")


@login_required
def upload_stl(request):
    if request.method == 'POST':
        print(request.FILES)

        upload_form = STLUploadForm(request.POST, request.FILES)
        process_form = STLProcessingForm(request.POST)

        if upload_form.is_valid() and process_form.is_valid():
            uploaded = upload_form.save(commit=False)
            user = request.user
            run_id = str(uuid.uuid4())
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            base_dir = os.path.join(settings.MEDIA_ROOT, 'stl_processed', f'{user.username}', f'run_{timestamp}_{run_id}')
            os.makedirs(base_dir, exist_ok=True)

            input_path = os.path.join(base_dir, "input.stl")
            output_filename = "output.ply"

            with open(input_path, 'wb+') as dest:
                for chunk in request.FILES['input_file'].chunks():
                    dest.write(chunk)

            try:
                simplify_face_count = process_form.cleaned_data.get('simplify_face_count')
                simplify_reduction_target = process_form.cleaned_data.get('simplify_reduction_target')

                metrics = process_stl(input_path, base_dir, simplify_face_count, simplify_reduction_target)

                for metric in metrics:
                    if isinstance(metrics[metric], np.ndarray):
                        metrics[metric] = metrics[metric].item()

                with open(os.path.join(base_dir, "metrics.json"), "w") as f:
                    json.dump(metrics, f)

                STLFile.objects.create(
                    user=user,
                    run_folder=base_dir,
                    input_file=input_path,
                    output_file=os.path.join(base_dir, output_filename),
                    metrics=metrics,
                    timestamp=timestamp
                )

                request.session['processed_stl_path'] = os.path.join(base_dir, output_filename)
                save_upload_log(user, request.FILES['input_file'].name, base_dir)

                return render(request, 'stl_viewer.html', {
                    'upload_form': upload_form,
                    'process_form': process_form,
                    'metrics': metrics,
                    'stl_filename': output_filename,
                    'upload_success': True
                })


            except Exception as e:
                return render(request, 'stl_viewer.html', {
                    'upload_form': upload_form,
                    'process_form': process_form,
                    'error': str(e),
                    'upload_success': False
                })
    else:
        upload_form = STLUploadForm()
        process_form = STLProcessingForm()

    return render(request, 'stl_viewer.html', {
        'upload_form': upload_form,
        'process_form': process_form,
        'upload_success': False
    })


@login_required
def recent_uploads(request):
    log_file = os.path.join(settings.MEDIA_ROOT, 'stl_processed', f'{request.user.username}', 'upload_log.txt')
    uploads = []

    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            for line in f:
                original, saved, time = line.strip().split(',')
                run_id = os.path.basename(saved)
                uploads.append({
                    'original': original,
                    'saved': saved,
                    'run_id': run_id,
                    'time': time,
                })
    # print(uploads)
    uploads.sort(key=lambda x: x['time'], reverse=True)

    return render(request, 'recent_uploads.html', {'uploads': uploads})


@login_required
def upload_details(request, run_id):
    user_dir = os.path.join(settings.MEDIA_ROOT, 'stl_processed', f'{request.user.username}')
    run_dir = os.path.join(user_dir, run_id)

    if not os.path.exists(run_dir):
        return HttpResponse("Run directory not found.", status=404)

    metrics_path = os.path.join(run_dir, "metrics.json")
    if os.path.exists(metrics_path):
        with open(metrics_path, "r") as f:
            metrics = json.load(f)
    else:
        metrics = {}

    return render(request, 'upload_details.html', {
        'metrics': metrics,
        'input_path': f"/media/stl_processed/{request.user.username}/{run_id}/input.stl",
        'output_path': f"/media/stl_processed/{request.user.username}/{run_id}/output.ply",
        'run_id': run_id,
    })


@login_required
def view_stl_3d(request):
    processed_path = request.session.get('processed_stl_path')
    if not processed_path or not os.path.exists(processed_path):
        return HttpResponse("Processed STL file not found. Please upload and process a file first.", status=404)

    try:
        import trimesh
        from trimesh.viewer import scene_to_html
        mesh = trimesh.load_mesh(processed_path)
        scene = trimesh.Scene(mesh)
        html_content = scene_to_html(scene)
        return HttpResponse(html_content)
    except Exception as e:
        return HttpResponse(f"Error rendering STL: {str(e)}", status=500)

from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.conf import settings
import os
import trimesh
from trimesh.viewer import scene_to_html

@login_required
def view_stl_3d_recent(request, run_id, file_type):
    """
    file_type: either 'input' or 'output'
    """
    user_id = request.user.username

    if file_type == 'input':
        filename = 'input.stl'
    elif file_type == 'output':
        filename = 'output.ply'
    else:
        return HttpResponse("Invalid file type requested.", status=400)

    file_path = os.path.join(
        settings.MEDIA_ROOT,
        'stl_processed',
        f'{user_id}',
        f'{run_id}',
        filename
    )

    if not os.path.exists(file_path):
        return HttpResponse(f"{filename} file not found for this run ID.", status=404)

    try:
        mesh = trimesh.load_mesh(file_path)
        scene = trimesh.Scene(mesh)
        html_content = scene_to_html(scene)
        return HttpResponse(html_content)
    except Exception as e:
        return HttpResponse(f"Error rendering file: {str(e)}", status=500)
