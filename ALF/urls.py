# Project's urls.py
from django.contrib import admin
from django.urls import path, include

urlpatterns = [
    path('admin/', admin.site.urls),
    path('home/', include('home.urls')),         # Home app
    path('dashboard/', include('dashboard.urls')),  # Dashboard app
    path('search/', include('search.urls')),  # Search app
]