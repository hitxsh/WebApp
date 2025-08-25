from django import forms
from django.contrib.auth.models import User
from django.contrib.auth.forms import UserCreationForm
from .models import UploadedImage

class SignUpForm(UserCreationForm):
    first_name = forms.CharField(max_length=30, required=True)
    last_name = forms.CharField(max_length=30, required=True)

    class Meta:
        model = User
        fields = ('username', 'first_name', 'last_name', 'password1', 'password2')

class ImageUploadForm(forms.ModelForm):
    class Meta:
        model = UploadedImage
        fields = ('image',)


# forms.py
from django import forms
from .models import STLFile

class STLUploadForm(forms.ModelForm):
    class Meta:
        model = STLFile
        fields = ['input_file']
        widgets = {
            'input_file': forms.FileInput(attrs={'accept': '.stl'})
        }

from django import forms
from django.core.exceptions import ValidationError

class STLProcessingForm(forms.Form):
    simplify_face_count = forms.IntegerField(
        required=False, 
        min_value=1, 
        label='Simplify Face Count',
        help_text='Target number of faces (leave empty if using reduction target)'
    )
    simplify_reduction_target = forms.FloatField(
        required=False, 
        min_value=0.0, 
        max_value=1.0, 
        label='Simplify Reduction Target (0-1)',
        help_text='Percentage to reduce mesh finally (leave empty if using face count)'
    )

    def clean(self):
        cleaned_data = super().clean()
        face_count = cleaned_data.get('simplify_face_count')
        reduction_target = cleaned_data.get('simplify_reduction_target')

        if face_count and reduction_target:
            raise ValidationError('Please provide only one of Simplify Face Count or Simplify Reduction Target, not both.')
        
        # Uncomment this if you want to make at least one field required
        # if not face_count and not reduction_target:
        #     raise ValidationError('Please provide either Simplify Face Count or Simplify Reduction Target.')

        return cleaned_data
