# Search's urls.py

from django.urls import path
from . import views

urlpatterns = [
    path('search/', views.search, name='search'),
    path('save_paper/', views.save_paper, name='save_paper'),
    path('papers/', views.papers, name='papers'),
    path('authors/', views.authors, name='authors'),
    path('home/', views.home, name='home'),
    path('dashboard/', views.dashboard, name='dashboard')
    # Add other search-related URLs here
]