from django.shortcuts import render
from django.db.models import Count, Sum
from search.models import Paper
import plotly.express as px
import pandas as pd
from .forms import PaperFilterForm

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

    return render(request, 'dashboard/dashboard.html', context)