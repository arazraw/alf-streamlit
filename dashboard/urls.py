# dashboard/urls.py
from django.urls import path
from . import views

app_name = 'dashboard'  # This is the application namespace

urlpatterns = [
    path('studies/', views.list_studies, name='list_studies'),
    # You can add more paths for additional views in the dashboard app here
]
