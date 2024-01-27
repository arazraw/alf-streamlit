# main urls.py

from django.urls import path
from . import views

urlpatterns = [
    path('search/', views.search, name='search'),
    path('save_paper/', views.save_paper, name='save_paper'),
    path('papers/', views.papers, name='papers'),
    path('authors/', views.authors, name='authors'),
    path('home/', views.home, name='home'),
    path('filter_papers/', views.filter_papers, name='filter_papers')
]
