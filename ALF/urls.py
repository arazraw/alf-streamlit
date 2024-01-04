# ALF/urls.py
from django.contrib import admin
from django.urls import include, path

urlpatterns = [
    path('admin/', admin.site.urls),
    path('dashboard/', include('dashboard.urls')),  # This includes the URLs from the dashboard app
    # Add more paths for other apps or functionalities here
]
