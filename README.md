Documentation for new template theme: https://github.com/phoenixcoded/nextro-free-bootstrap-admin-template

PubMex XML file format info: https://www.ncbi.nlm.nih.gov/books/NBK3828/#publisherhelp.PubMed_XML_Tagged_Format   

Semantic scholar API: https://api.semanticscholar.org/api-docs/#tag/Author-Data/operation/get_graph_get_author_papers  

Now the view "Papers" will be updated every time we run the command "python manage.py update_papers.py". This command script will loop through all doi in Main and check if they exist in PApers, otherwise it will fetch the records from PubMed using the functions created in /search/views.py. We should be able to do this also for the Authors table.
