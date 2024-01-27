# dash_apps.py inside the vetu app

import dash
from dash import html, dcc
import pandas as pd
import plotly.express as px
from django_plotly_dash import DjangoDash
from .models import Paper
import plotly.graph_objects as go

def create_papers_dash_app():
    # Create a Dash application
    app = DjangoDash('PapersVisualization')

    # Fetch data from your Paper model
    papers_data = Paper.objects.all().values('akademisk')
    df = pd.DataFrame(papers_data)

    # Count the number of papers for each value of 'akademisk'
    akademisk_counts = df['akademisk'].value_counts()

    # Create a pie chart to visualize the counts
    fig = go.Figure(data=[go.Pie(labels=akademisk_counts.index, values=akademisk_counts.values)])
    fig.update_layout(
        autosize=True,
        margin=dict(l=20, r=20, t=40, b=20),
    )

    graph_config = {'displayModeBar': True, 'responsive': True}

    # Set up the layout of the Dash app
    app.layout = html.Div([
    html.H1('Graph title'),
    dcc.Graph(
        figure=fig,
        config=graph_config,
        style={'height': '100vh'}  # 100% of the viewport height
    ),
])

    return app
