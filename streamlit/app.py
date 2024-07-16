import streamlit as st
import pandas as pd
import psycopg2
import matplotlib.pyplot as plt
import plotly.express as px
from supabase import create_client, Client

file_path_university = '/Users/xanerc/Documents/Vetu/alf/RegionerAkademier/affiliations_university_norm.csv'
universities = pd.read_csv(file_path_university)
universities['Code'] = universities['Code'].astype(str)

file_path_university2 = '/Users/xanerc/Downloads/affiliations_university_decoder_list.csv'
universities2 = pd.read_csv(file_path_university2)
universities2['Code'] = universities2['Code'].astype(str)


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

def plot_line_graph(df, x_column, y_column, title):
    fig = px.line(df, x=x_column, y=y_column, title=title)
    fig.update_layout(title=title)
    return fig

# Set the page to wide mode
st.set_page_config(layout="wide")

# Inject custom CSS to style the radio buttons as links and increase text size
st.markdown("""
    <style>
    div.stRadio > div {
        display: flex;
        flex-direction: column;
    }
    div.stRadio > div > label {
        display: block;
        padding: 10px 0;
        cursor: pointer;
        text-decoration: none;
        color: #007bff;
        font-weight: bold;
        font-size: 42px;
    }
    div.stRadio > div > label > input {
        display: none;
    }
    div.stRadio > div > label > div:nth-of-type(1) {
        display: none;
    }
    </style>
    """, unsafe_allow_html=True)

# Create a sidebar
st.sidebar.header('Välj verktyg')

# Navigation menu
navigation = st.sidebar.radio('', ('Hem', 'Översikt', 'Akademi & Högskola', 'Region (ALF)', 'Tidsskrifter', 'Forskare', 'Finansiärer', 'Innovation', 'Sök Artiklar'))

st.title(navigation)

if navigation == 'Hem':
    st.subheader("Welcome!")

elif navigation == 'Översikt':
    
    col1, col2 = st.columns(2)

    with col1:
        # Function to fetch data from the database
        def fetch_impact_data():
            conn = create_conn()
            query = "SELECT citations, impactful_citations FROM vetu_impact"
            df = pd.read_sql_query(query, conn)
            conn.close()
            return df

        # Fetch data
        impact_df = fetch_impact_data()

        # Calculate totals
        total_citations = impact_df['citations'].sum()
        total_impactful_citations = impact_df['impactful_citations'].sum()
        total_papers = impact_df['citations'].count()

        # Display the info box
        st.markdown(f"""
        <div style="background-color: #f9f9f9; padding: 20px; border-radius: 5px;">
            <h3 style="color: #333;">Impact Summary</h3>
            <p><strong>Total Citations:</strong> {total_citations}</p>
            <p><strong>Total Impactful Citations:</strong> {total_impactful_citations}</p>
            <p><strong>Total Number of Papers:</strong> {total_papers}</p>
        </div>
        """, unsafe_allow_html=True)

        # Add a small space
        st.write('---')

    with col2:
        # Function to fetch data from the database
        def fetch_papers_per_year():
            conn = create_conn()
            query = "SELECT year FROM vetu_paper"
            df = pd.read_sql_query(query, conn)
            conn.close()
            return df

        # Fetch data
        df = fetch_papers_per_year()

        # Group by year and count the number of papers
        papers_per_year = df['year'].value_counts().reset_index()
        papers_per_year.columns = ['Year', 'Total Papers']
        papers_per_year = papers_per_year.sort_values('Year')

        # Plot the bar chart
        fig = px.bar(papers_per_year, x='Year', y='Total Papers', title='Total Papers Published Each Year',
                    labels={'Year': 'Year', 'Total Papers': 'Number of Papers'}
                    )
        
        # Update the x-axis range to start at 1990
        fig.update_layout(
            xaxis=dict(
                range=[1989, 2025],
                tickmode='linear',
                tick0=1990,
                dtick=5
            )
        )

        # Display the plot in Streamlit
        st.plotly_chart(fig)

    # Function to fetch data from the database
    def fetch_citation_data():
        conn = create_conn()
        query = "SELECT year, citations FROM vetu_impact"
        df = pd.read_sql_query(query, conn)
        conn.close()
        return df

    # Fetch data
    citation_df = fetch_citation_data()

    # Categorize papers based on citation counts
    citation_df['citation_category'] = pd.cut(
        citation_df['citations'],
        bins=[-1, 1, 5, citation_df['citations'].max()],
        labels=['≤1 Citation', '2-5 Citations', '>5 Citations']
    )

    # Calculate the percentage composition for each year
    composition_df = citation_df.groupby(['year', 'citation_category']).size().reset_index(name='count')
    total_per_year = composition_df.groupby('year')['count'].transform('sum')
    composition_df['percentage'] = composition_df['count'] / total_per_year * 100

    # Pivot the data for plotting
    pivot_df = composition_df.pivot(index='year', columns='citation_category', values='percentage').fillna(0)

    # Plot the percentage composition as a stacked bar chart
    fig = px.bar(pivot_df, x=pivot_df.index, y=pivot_df.columns, title='Percentage Composition of Papers by Citation Counts (1990-2024)',
                labels={'value': 'Percentage', 'year': 'Year'}, 
                barmode='stack')

    # Update the x-axis range to start at 1990
    fig.update_layout(
        xaxis=dict(
            range=[1990, 2024],
            tickmode='linear',
            tick0=1990,
            dtick=1
        ),
        yaxis=dict(
            title='Percentage'
        )
    )

    # Display the plot in Streamlit
    st.plotly_chart(fig)

elif navigation == 'Akademi & Högskola':

    # Function to fetch data from the database based on filters
    def fetch_data(university, institute, department, from_year, to_year):
        conditions = []
        if university != "All":
            conditions.append(f"affiliations LIKE '%{university}%'")
        if institute != "All":
            conditions.append(f"affiliations LIKE '%{institute}%'")
        if department != "All":
            conditions.append(f"affiliations LIKE '%{department}%'")
        conditions.append(f"year >= {from_year}")
        conditions.append(f"year <= {to_year}")

        where_clause = " AND ".join(conditions)

        query = f"""
            SELECT year, COUNT(*) as publication_count
            FROM vetu_paper
            WHERE {where_clause}
            GROUP BY year
            ORDER BY year;
        """
        conn = create_conn()
        df = get_data(conn, query)
        conn.close()
        return df

    # Create a box containing four dropdown menus
    with st.container():
        # Generate a list of years from 1990 to 2024
        fran_ar_list = list(range(1990, 2025))
        fran_ar = 1990

        # Arrange dropdown menus in columns
        col1, col2 = st.columns(2)
        col3, col4, col5 = st.columns(3)
        col6, col7 = st.columns(2)
        col8, col9, col10 = st.columns(3)

        # Other dropdown menus
        with col1:
            year_range = st.slider('År:', min_value=1990, max_value=2024, value=(1990, 2024)) # År slider
            fran_ar, till_ar = year_range

        with col2:
            pass

        with col3:
            selected_university = st.selectbox('Universitet:', ["All"] + universities2[universities2['Code'].str.count('\.') == 0]['Department'].tolist(), index=0) # Universitet
            if selected_university != "All":
                selected_university_code = universities2[universities2['Department'] == selected_university]['Code'].values[0]
            else:
                selected_university_code = ""

        with col4:
            if selected_university == "ALL":
                st.selectbox('Institut:', ["All"])
            else:
                selected_institute = st.selectbox('Institut:',
                ["All"] + universities2[
                    (universities2['Code'].str.startswith(selected_university_code + '.')) & (universities2['Code'].str.count('\.')== 1)]['Department'].tolist(), index=0
                ) # Institut
                if selected_institute != "All":
                    selected_institute_code = universities2[universities2['Department'] == selected_institute]['Code'].values[0]
                else:
                    selected_institute_code = ""

        with col5:    
            if selected_institute == "ALL":
                st.selectbox('Department:', ["All"])
            else:
                selected_department = st.selectbox('Department:', 
                ["All"] + universities2[
                    (universities2['Code'].str.startswith(selected_institute_code + '.')) & (universities2['Code'].str.count('\.') == 2)]['Department'].tolist(), index=0
                ) # Avdelning

        #with col5:
            #selected_topic = st.selectbox('Ämne:', ['Option 2 - A', 'Option 2 - B', 'Option 2 - C']) # Ämne

        with col6:
            if st.checkbox('Compare'):
                with col8:
                    selected_university_comp = st.selectbox('Jämför Universitet:', ["All"] + universities2[universities2['Code'].str.count('\.') == 0]['Department'].tolist(), index=0) # Universitet
                    if selected_university_comp != "All":
                        selected_university_code_comp = universities2[universities2['Department'] == selected_university_comp]['Code'].values[0]
                    else:
                        selected_university_code_comp = ""

                with col9:
                    if selected_university_comp == "ALL":
                        st.selectbox('Institut:', ["All"])
                    else:
                        selected_institute_comp = st.selectbox('Jämför Institut:',
                        ["All"] + universities2[
                            (universities2['Code'].str.startswith(selected_university_code_comp + '.')) & (universities2['Code'].str.count('\.')== 1)]['Department'].tolist(), index=0
                        ) # Institut
                        if selected_institute_comp != "All":
                            selected_institute_code_comp = universities2[universities2['Department'] == selected_institute_comp]['Code'].values[0]
                        else:
                            selected_institute_code_comp = ""

                with col10:    
                    if selected_institute_comp == "ALL":
                        st.selectbox('Department:', ["All"])
                    else:
                        selected_department_comp = st.selectbox('Jämför Department:', 
                        ["All"] + universities2[
                            (universities2['Code'].str.startswith(selected_institute_code_comp + '.')) & (universities2['Code'].str.count('\.') == 2)]['Department'].tolist(), index=0
                        ) # Avdelning
                
                data2 = fetch_data(selected_university_comp, selected_institute_comp, selected_department_comp, fran_ar, till_ar)

            else:
                data2 = pd.DataFrame()

        with col7:
            pass



    # Function to fetch data from the database based on filters
        def fetch_data(university, institute, department, from_year, to_year):
            conditions = []
            if university != "All":
                conditions.append(f"affiliations LIKE '%{university}%'")
            if institute != "All":
                conditions.append(f"affiliations LIKE '%{institute}%'")
            if department != "All":
                conditions.append(f"affiliations LIKE '%{department}%'")
            conditions.append(f"year >= {from_year}")
            conditions.append(f"year <= {to_year}")

            where_clause = " AND ".join(conditions)

            query = f"""
                SELECT year, COUNT(*) as publication_count
                FROM vetu_paper
                WHERE {where_clause}
                GROUP BY year
                ORDER BY year;
            """
            conn = create_conn()
            df = get_data(conn, query)
            conn.close()
            return df

    # Fetch the data
    data = fetch_data(selected_university, selected_institute, selected_department, fran_ar, till_ar)

    def plot_data(data):
        if data.empty:
            st.write("No data available for the selected criteria.")
        else:
            fig = px.bar(data, x='year', y='publication_count', title='Publications Over Time',
              labels={'year': 'Year', 'publication_count': 'Number of Publications'})
            fig.update_layout(
            xaxis=dict(
                tickmode='linear',
                tick0=fran_ar,
                dtick=1,
                range=[fran_ar-0.5, till_ar+0.5])  # Use selected from_year and to_year for range
            )
            st.plotly_chart(fig)

    # Plot the data
    if data2.empty:
        # Plot only data
        fig = px.bar(data, x='year', y='publication_count', title='Publications Over Time',
                    labels={'year': 'Year', 'publication_count': 'Number of Publications'})
        fig.update_layout(
            xaxis=dict(
                tickmode='linear',
                tick0=data['year'].min(),
                dtick=1
            )
        )
        st.plotly_chart(fig)

    elif not data.empty and not data2.empty:
        # Combine data for side-by-side plotting
        data['Type'] = f"{selected_university} - {selected_institute} - {selected_department}"
        data2['Type'] = f"{selected_university_comp} - {selected_institute_comp} - {selected_department_comp}"

        combined_data = pd.concat([data, data2])

        fig = px.bar(combined_data, x='year', y='publication_count', color='Type', barmode='group',
                    title='Publications Over Time',
                    labels={'year': 'Year', 'publication_count': 'Number of Publications'})
        fig.update_layout(
            xaxis=dict(
                tickmode='linear',
                tick0=min(data['year'].min(), data2['year'].min()),
                dtick=1
            ),
            legend=dict(
            orientation='h',  # Horizontal legend
            yanchor='top',  # Anchor the legend at the top
            y=-0.2,  # Position the legend below the graph
            xanchor='center',  # Center the legend horizontally
            x=0.5)  # Align the legend at the center of the x-axis
        )
        st.plotly_chart(fig)


elif navigation == 'Finansiärer':
    st.title('Överskådliggöra')

elif navigation == 'Tidsskrifter':

    # Create a connection to the database
    conn = create_conn()

    # SQL query to get unique affiliations
    affiliations_query = "SELECT DISTINCT affiliations FROM vetu_paper;"

    # Get unique affiliations
    affiliations_df = get_data(conn, affiliations_query)
    affiliations = affiliations_df['affiliations'].tolist()

    # Create a dropdown menu for affiliations
    selected_affiliation = st.selectbox('Select an Affiliation:', affiliations)

    # SQL query to get the count of publications per year for the selected affiliation
    query = f"""
        SELECT year, COUNT(*) AS publication_count
        FROM vetu_paper
        WHERE affiliations = '{selected_affiliation}'
        GROUP BY year
        ORDER BY year;
    """

    # Get data based on the selected affiliation
    df_year = get_data(conn, query)

    # Plot the data as a line graph
    st.plotly_chart(plot_line_graph(df_year, 'year', 'publication_count', f'Number of Publications per Year for {selected_affiliation}'))

elif navigation == 'Forskare':

    # Function to fetch data from the database
    def fetch_data_forskare():
        conn = create_conn()
        query = "SELECT * FROM vetu_authorimpact"
        df = pd.read_sql_query(query, conn)
        conn.close()
        return df

    # Fetch data
    df = fetch_data_forskare()

    # Create search bar
    search_query = st.text_input("Search within data", "")

    # Filter the DataFrame based on the search query
    if search_query:
        filtered_df = df[df.apply(lambda row: row.astype(str).str.contains(search_query, case=False).any(), axis=1)]
    else:
        filtered_df = df

    # Rename the columns
    filtered_df = filtered_df.rename(columns={
        "name": "Author",
        "citations": "Total Citations",
        "impactful_citations": "Impactful Citations",
        "paper_count": "Number of Papers",
        "affiliations": "Affiliation"
    })
    # Add a unique identifier for each row
    filtered_df["Unique Author"] = filtered_df["Author"] + " (" + filtered_df.index.astype(str) + ")"

    # Reset the index and drop the original index
    filtered_df = filtered_df.sort_values(by="Total Citations", ascending=False)
    filtered_df = filtered_df.reset_index(drop=True)

    # Select relevant columns for display
    columns_to_display = ["Author", "Total Citations", "Impactful Citations", "Number of Papers", "Affiliation"]
    filtered_df1 = filtered_df[columns_to_display]

    # Display filtered DataFrame
    st.dataframe(filtered_df1)

    # Sort by "Total Citations" in descending order and select the top 10
    filtered_df2 = filtered_df.sort_values(by="Total Citations", ascending=False).head(20)

    # Example plot (optional)
    if not filtered_df2.empty:
        fig = px.bar(filtered_df2, x='Unique Author', y='Total Citations', title='Citations by Author',
                    labels={'Unique Author': 'Author', 'Citations': 'Number of Citations'})
        st.plotly_chart(fig)
  

elif navigation == 'Innovation':
    st.title('Överskådliggöra')

elif navigation == 'Region (ALF)':
    st.title('Överskådliggöra')

elif navigation == 'Sök Artiklar':
    st.title('Överskådliggöra')

