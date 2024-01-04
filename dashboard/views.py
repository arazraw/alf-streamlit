# Dashboard's views.py
from django.shortcuts import render
from django.db import models
from .models import Main
import plotly.express as px
import pandas as pd

def dashboard(request):
    # Query the Main model to get the number of papers per year
    data = Main.objects.values('year').annotate(paper_count=models.Count('year')).order_by('year')

    # Create a DataFrame from the queryset
    df = pd.DataFrame.from_records(data)

    # Create a Plotly bar chart
    fig = px.bar(df, x='year', y='paper_count', labels={'paper_count': 'Number of Papers'}, title='Number of Papers per Year')
    
    # Convert the Plotly chart to HTML
    chart_div = fig.to_html(full_html=False, default_height=500, default_width=800)

    return render(request, 'dashboard/dashboard.html', {'chart_div': chart_div})
