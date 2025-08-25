import random
from PIL import Image, ImageDraw

# Image dimensions
width, height = 500, 500

# Create a new black background image
image = Image.new("RGB", (width, height), "black")
draw = ImageDraw.Draw(image)

# Define the number of points for the polygon (minimum 3 points)
num_points = random.randint(3, 10)

# Generate random points within the image dimensions
points = [(random.randint(0, width), random.randint(0, height)) for _ in range(num_points)]

# Draw the polygon in yellow
draw.polygon(points, fill="yellow")

# Save the image as a TIFF file
image.save("yellow_polygon1.tif")

print("Image saved as yellow_polygon.tif")
