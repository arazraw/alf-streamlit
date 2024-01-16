# Search's views.py
from django.shortcuts import render
from vetu.models import Paper, Author, Impact  # Update this line
from django.views.decorators.http import require_POST
from django.http import JsonResponse
from django.db.models import Count, Sum
import plotly.express as px
import pandas as pd
from .forms import PaperFilterForm

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

        redirect_url = request.META.get('HTTP_REFERER', '/')
        print(redirect_url)
        response_data = {'status': 'success', 'redirect_url': redirect_url}
        return JsonResponse(response_data)

    return JsonResponse({'status': 'error form type error'}, status=400)

def dashboard(request):
    # Handle the form submission
    if request.method == 'POST':
        form = PaperFilterForm(request.POST)
        if form.is_valid():
            start_year = form.cleaned_data['start_year']
            end_year = form.cleaned_data['end_year']
            if start_year and end_year:
                papers = Paper.objects.filter(year__range=(start_year, end_year))
            else:
                papers = Paper.objects.all()
    else:
        form = PaperFilterForm()
        papers = Paper.objects.all()

    # Calculate the total number of papers and total citations
    paper_count = papers.count()
    total_citations = papers.aggregate(total_citations=Sum('citations'))['total_citations']

    # Query for the number of papers per year
    papers_per_year_data = papers.values('year').annotate(paper_count=Count('year')).order_by('year')
    # Query for the total number of citations per year
    citations_per_year_data = papers.values('year').annotate(total_citations=Sum('citations')).order_by('year')

    # Papers per Year Chart
    if papers_per_year_data:
        df_papers = pd.DataFrame.from_records(papers_per_year_data)
        fig_papers = px.bar(df_papers, x='year', y='paper_count', labels={'paper_count': 'Number of Papers'}, title='Number of Papers per Year')
        chart_papers_div = fig_papers.to_html(full_html=False, default_height=500, default_width=800)
    else:
        chart_papers_div = 'No data available for papers per year'

    # Citations per Year Chart
    if citations_per_year_data:
        df_citations = pd.DataFrame.from_records(citations_per_year_data)
        fig_citations = px.bar(df_citations, x='year', y='total_citations', labels={'total_citations': 'Total Citations'}, title='Total Citations per Year')
        chart_citations_div = fig_citations.to_html(full_html=False, default_height=500, default_width=800)
    else:
        chart_citations_div = 'No data available for citations per year'

    context = {
        'form': form,
        'chart_papers_div': chart_papers_div,
        'chart_citations_div': chart_citations_div,
        'paper_count': paper_count,
        'total_citations': total_citations
    }

    return render(request, 'vetu/dashboard.html', context)


def papers(request):
    papers = Paper.objects.all()
    return render(request, 'vetu/papers.html', {'papers': papers})

def authors(request):
    authors = Author.objects.all()
    return render(request, 'vetu/authors.html', {'authors': authors})

def home(request):
    return render(request, 'vetu/home.html')
