Updates:
1. New template used for frontend: Nextro Able. Documentation: https://github.com/phoenixcoded/nextro-free-bootstrap-admin-template
2. Added folder "SwePub" with all medical reserach articles from 1980 to 2024 in SwePub, and a script to parse JSON to a dataframe. These papers (>400k) can be our initial starting point for the app when it goes live.

PubMex XML file format info: https://www.ncbi.nlm.nih.gov/books/NBK3828/#publisherhelp.PubMed_XML_Tagged_Format   

Semantic scholar API: https://api.semanticscholar.org/api-docs/#tag/Author-Data/operation/get_graph_get_author_papers  

Now the view "Papers" will be updated every time we run the command "python manage.py update_papers.py". This command script will loop through all doi in Main and check if they exist in PApers, otherwise it will fetch the records from PubMed using the functions created in /search/views.py. We should be able to do this also for the Authors table.
