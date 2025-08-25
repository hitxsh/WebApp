from django.urls import path
from . import views
from django.conf import settings
from django.conf.urls.static import static

urlpatterns = [
    path('', views.home, name='home'),
    path('login/', views.login_view, name='login'),
    path('logout/', views.logout_view, name='logout'),
    path('signup/', views.signup_view, name='signup'),
    path('upload/', views.upload_image, name='upload'),
    path('recents/', views.recent_results, name='recents'),
    # path('upload_stl/', views.upload_stl, name='upload_stl'),
    # path('view_stl/', views.view_stl, name='view_stl'),
    # path('view_stl/<path:stl_path>/', views.view_stl, name='view_stl_with_path'),
    path('upload_stl/', views.upload_stl, name='upload_stl'),
    path('view_stl_3d/', views.view_stl_3d, name='view_stl_3d'),
    path('recent/', views.recent_uploads, name='recent_uploads'),
    path('recent/<str:run_id>/', views.upload_details, name='upload_details'),
    path('view_stl_3d_recent/<str:run_id>/<str:file_type>/', views.view_stl_3d_recent, name='view_stl_3d_recent'),

]+ static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)


if settings.DEBUG:
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
