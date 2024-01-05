# Search's views.py
from django.shortcuts import render
from search.models import Paper, Author  # Update this line
from django.views.decorators.http import require_POST
from django.http import JsonResponse

# Importing scripts to fetch and save articles
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
        print(data)  # Add this line
        # Check if each article is already saved
        for record in data:
            record['saved'] = Paper.objects.filter(doi=record['DOI']).exists()

    return render(request, 'search/search.html', {'data': data})

@require_POST
def save_paper(request):
    if request.POST.get('form_type') == 'save_paper':
        title = request.POST.get('title')
        doi = request.POST.get('doi')
        year = request.POST.get('year')
        pmid = request.POST.get('pmid')
        abstract_text = request.POST.get('abstract_text')
        journal_title = request.POST.get('journal_title')
        publication_type = request.POST.get('publication_type')
        month = request.POST.get('month')
        
        # Validate the 'year' field
        try:
            year = int(year)  # Try to convert 'year' to an integer
        except ValueError:
            # Handle the case where 'year' is not a valid integer
            year = None  # Set a default value or None depending on your needs

        # Create a new Paper object and save it to the database
        new_paper = Paper(
            title=title, 
            doi=doi, 
            year=year, 
            pmid=pmid,
            abstract_text=abstract_text,
            journal_title=journal_title,
            publication_type=publication_type,
            month=month,
        )
        new_paper.save()
        
        # Return a JSON response
        return JsonResponse({'status': 'success'})

    return JsonResponse({'status': 'error'}, status=400)

def papers(request):
    papers = Paper.objects.all()
    return render(request, 'search/papers.html', {'papers': papers})

def authors(request):
    authors = Author.objects.all()
    return render(request, 'search/authors.html', {'authors': authors})