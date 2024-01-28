import streamlit as st
import pandas as pd
import psycopg2
import matplotlib.pyplot as plt
import plotly.express as px

# Set the page to wide mode
st.set_page_config(layout="wide")

def create_conn():
    conn = psycopg2.connect(
        host='aws-0-eu-central-1.pooler.supabase.com',
        dbname='postgres',
        user='postgres.tszoiepgwjrdmbviywxz',
        password='wBv(2+ayMwAT5_N',
        port='6543'
    )
    return conn

def get_data(conn, query):
    df = pd.read_sql(query, conn)
    return df

def plot_distribution(df, column, title):
    fig = px.histogram(df, x=column)
    fig.update_layout(title=title)
    return fig

def main():
    st.title("Vetenskaplig utvärdering av svensk medicinsk forskning")

    conn = create_conn()

    # Example plot 1: Distribution of publication years
    query = "SELECT year FROM vetu_paper;"
    df_year = get_data(conn, query)
    st.plotly_chart(plot_distribution(df_year, 'year', 'Distribution of Publication Years'))

    # Add more plots here using different queries and columns from your database
    # Example: 
    # query = "SELECT topic FROM vetu_paper;"
    # df_topic = get_data(conn, query)
    # st.plotly_chart(plot_distribution(df_topic, 'topic', 'Distribution of Topics'))

    # More plots can be added similarly

if __name__ == "__main__":
    main()
