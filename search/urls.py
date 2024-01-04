# Search's urls.py

from django.urls import path
from . import views

urlpatterns = [
    path('', views.search, name='search'),
    # Add other search-related URLs here
]
