# Search's views.py
from django.shortcuts import render
from vetu.models import Paper, Author
from django.views.decorators.http import require_POST
from django.http import JsonResponse
from django.db.models import Count, Sum
import plotly.express as px
import pandas as pd
from .forms import PaperFilterForm
from django.core import serializers

# Importing scripts to fetch and save articles
from .pubmed_script import fetch_articles

def papers(request):
    papers = Paper.objects.all()
    return render(request, 'vetu/papers.html', {'papers': papers})

def authors(request):
    authors = Author.objects.all()
    return render(request, 'vetu/authors.html', {'authors': authors})

def home(request):
    return render(request, 'vetu/home.html')


def search(request):
    data = []
    if request.method == 'POST' and request.POST.get('form_type') == 'search':
        term = request.POST.get('term')
        total_articles = int(request.POST.get('total_articles'))
        mindate = request.POST.get('mindate')
        maxdate = request.POST.get('maxdate')

        # Fetch articles using your script
        data = fetch_articles(term, total_articles, mindate, maxdate)
        # print(data)  # Add this line
        # Check if each article is already saved
        for record in data:
            record['saved'] = Paper.objects.filter(doi=record['DOI']).exists()

    return render(request, 'vetu/search.html', {'data': data})

@require_POST
def save_paper(request):
    # print(request.POST.get('form_type'))
    # print(request.POST.get('title'))
    # print(request.POST.get('doi'))
    # print(request.POST.get('year'))
    # print(request.POST.get('pmid'))
    # print(request.POST.get('abstract_text'))
    # print(request.POST.get('journal_title'))
    # print(request.POST.get('publication_type'))
    # print(request.POST.get('month'))
    if request.POST.get('form_type') == 'save_paper':
        title = request.POST.get('title')
        doi = request.POST.get('doi')
        year = request.POST.get('year')
        pmid = request.POST.get('pmid')
        abstract_text = request.POST.get('abstract_text')
        journal_title = request.POST.get('journal_title')
        publication_type = request.POST.get('publication_type')
        month = request.POST.get('month')
        affiliations = request.POST.get('affiliations')
        
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
            affiliations=affiliations,
        )
        new_paper.save()

        redirect_url = request.META.get('HTTP_REFERER', '/')
        print(redirect_url)
        response_data = {'status': 'success', 'redirect_url': redirect_url}
        return JsonResponse(response_data)

    return JsonResponse({'status': 'error form type error'}, status=400)


def filter_papers(request):
    print("filter_papers view called")  # Add this line
    start_year = request.GET.get('start_year')
    end_year = request.GET.get('end_year')
    akademisk = request.GET.get('akademisk') == 'on'
    region = request.GET.get('region') == 'on'

    papers = Paper.objects.all()

    if start_year and end_year:
        start_year = int(start_year)
        end_year = int(end_year)
        papers = papers.filter(year__range=(start_year, end_year))
    if akademisk:
        papers = papers.filter(akademisk=akademisk)
    if region:
        papers = papers.filter(region=region)

    # Check if any filters are applied
    if not (start_year or end_year or akademisk or region):
        # No filters applied, return first 20 papers
        papers = papers[:20]

    return render(request, 'vetu/papers_table.html', {'papers': papers})

