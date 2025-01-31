## Vår skiss på layout till översikten
- https://www.figma.com/file/F5LnbHafAd6FdruQYRXWzb/FOUI?type=whiteboard&node-id=0-1&t=cCrRPL2h3OGHPwtd-0

# Important / useful resources
- Härifrån hämtar vi merparten av all information om artiklarna: https://api.semanticscholar.org/api-docs/graph#tag/Paper-Data/operation/post_graph_get_papers
- Collection of resources: https://github.com/topics/citation-analysis
- Seems important: https://github.com/scholarly-python-package/scholarly
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10132428/
- PubMex XML file format info: https://www.ncbi.nlm.nih.gov/books/NBK3828/#publisherhelp.PubMed_XML_Tagged_Format   
- Semantic scholar API: https://api.semanticscholar.org/api-docs/#tag/Author-Data/operation/get_graph_get_author_papers  

## Starting with a bulk of paper
Added folder "SwePub" with all medical reserach articles from 1980 to 2024 in SwePub, and a script to parse JSON to a dataframe. These papers (>400k) can be our initial starting point for the app when it goes live.

Currently the command import_swepub has a hardcoded destination for saving the data, this is to avoid saving it into the project as the file is 37GB when uncompacted. When attempting to run ensure that this is individualy set to desired destination. There exists code to later define a path within the project.

## Commands we created
- The view "Papers" will be updated every time we run the command "python manage.py update_papers.py". This command script will loop through all doi in Main and check if they exist in Papers, otherwise it will fetch the records from PubMed using the functions created in /search/views.py.
- There are now multiple commands, there are those to update_authors, update_impact, fill_affiliations also does binary checks for certain types of affiliations.
- In pipeline(see dev branch) there are several more which get, proccess and load data from swepub, PubMed and semantic schloar

## Affiliations
- The affiliation codes found in the papers column correspond to universities and departments in sweden, the decoding can be found in RegionerAkademier folder, it is currently saved in code form for easier handling.
