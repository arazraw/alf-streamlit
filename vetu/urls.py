# vetu app's urls.py

from django.urls import path
from . import views

urlpatterns = [
    path('search/', views.search, name='search'),
    path('save_paper/', views.save_paper, name='save_paper'),
    path('papers/', views.papers, name='papers'),
    path('authors/', views.authors, name='authors'),
    path('home/', views.home, name='home'),
    path('filter_papers/', views.filter_papers, name='filter_papers'),
    path('papers_dash/', views.papers_dash_view, name='papers_dash'),
    # Remove the django_plotly_dash line from here
]
