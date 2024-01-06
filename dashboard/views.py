from django.shortcuts import render
from django.db.models import Count, Sum
from search.models import Paper
import plotly.express as px
import pandas as pd

def dashboard(request):
    # Query for the number of papers per year
    papers_per_year_data = Paper.objects.values('year').annotate(paper_count=Count('year')).order_by('year')
    # Query for the total number of citations per year
    citations_per_year_data = Paper.objects.values('year').annotate(total_citations=Sum('citations')).order_by('year')

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
        'chart_papers_div': chart_papers_div,
        'chart_citations_div': chart_citations_div
    }

    return render(request, 'dashboard/dashboard.html', context)
