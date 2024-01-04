from django.shortcuts import render

# Create your views here.
from django.shortcuts import render
from .models import Main

def list_studies(request):
    studies = Main.objects.all()  # Retrieve all studies
    return render(request, 'dashboard/list_studies.html', {'studies': studies})
