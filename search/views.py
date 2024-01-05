# Search's views.py
from django.shortcuts import render
from django.shortcuts import redirect
#from django.views.decorators.csrf import csrf_exempt
from dashboard.models import Main
from django.views.decorators.http import require_POST
from .pubmed_script import fetch_articles

def search(request):
    data = []
    if request.method == 'POST' and request.POST.get('form_type') == 'search':
        term = request.POST.get('term')
        total_articles = int(request.POST.get('total_articles'))
        mindate = request.POST.get('mindate')
        maxdate = request.POST.get('maxdate')

        # Fetch articles using your script
        data = fetch_articles(term, total_articles, mindate, maxdate)

    return render(request, 'search/search.html', {'data': data})

from django.contrib import messages



from django.http import JsonResponse
from django.views.decorators.http import require_POST

@require_POST
def save_paper(request):
    if request.POST.get('form_type') == 'save_paper':
        title = request.POST.get('title')
        doi = request.POST.get('doi')
        year = request.POST.get('year')
        pmid = request.POST.get('pmid')
        
        # Validate the 'year' field
        try:
            year = int(year)  # Try to convert 'year' to an integer
        except ValueError:
            # Handle the case where 'year' is not a valid integer
            year = None  # Set a default value or None depending on your needs

        # Create a new Main object and save it to the database
        new_main = Main(title=title, doi=doi, year=year, pmid=pmid)
        new_main.save()
        
        # Return a JSON response
        return JsonResponse({'status': 'success'})

    return JsonResponse({'status': 'error'}, status=400)
