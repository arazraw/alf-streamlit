# Dashboard's urls.py
from django.urls import path
from . import views

urlpatterns = [
    path('', views.dashboard, name='dashboard'),  # Main dashboard view
    # Add other dashboard-related URLs here
]
