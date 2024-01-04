# Search's urls.py

from django.urls import path
from . import views

urlpatterns = [
    path('', views.search, name='search'),
    path('save_paper/', views.save_paper, name='save_paper'),
    # Add other search-related URLs here
]
