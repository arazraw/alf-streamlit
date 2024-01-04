# Dashboard's views.py

from django.shortcuts import render

def dashboard(request):
    # Your dashboard view logic goes here
    return render(request, 'dashboard/dashboard.html')
