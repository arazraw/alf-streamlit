# Dashboard's urls.py
from django.urls import path
from . import views

urlpatterns = [
    path('', views.dashboard, name='dashboard'),
    # Add other dashboard-related URLs here
]
