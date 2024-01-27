# Project's urls.py
from django.contrib import admin
from django.urls import path, include
from django.views.generic import RedirectView

urlpatterns = [
    path('', RedirectView.as_view(url='/home/', permanent=True)),  # Redirect from root to /home/
    path('admin/', admin.site.urls),
    path('', include('vetu.urls')),  # This includes 'vetu.urls' at the root path - could conflict with the redirect above
    path('django_plotly_dash/', include('django_plotly_dash.urls')),  # Include this line for Django Plotly Dash
]