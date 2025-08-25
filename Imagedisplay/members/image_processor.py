from PIL import Image
import os
from django.conf import settings

def get_image_dimensions(image_path):
    """
    Takes an image file path, opens it, and returns its dimensions (width, height).
    """
    full_path = os.path.join(settings.MEDIA_ROOT, image_path)

    try:
        with Image.open(full_path) as img:
            width, height = img.size
        return width, height
    except Exception as e:
        return f"Error processing image: {e}"
