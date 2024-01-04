# Search's urls.py

from django.shortcuts import render

#def search(request):
#    return render(request, 'search/search.html')


# Import your PubMed script here
from .pubmed_script import fetch_articles

def search(request):
    data = []
    if request.method == 'POST':
        term = request.POST.get('term')
        total_articles = int(request.POST.get('total_articles'))
        mindate = request.POST.get('mindate')
        maxdate = request.POST.get('maxdate')

        # Fetch articles using your script
        data = fetch_articles(term, total_articles, mindate, maxdate)

    return render(request, 'search/search.html', {'data': data})
