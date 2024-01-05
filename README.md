Main search view must enable admins to search for articles in PubMed (and later also preprintservers, Web of Science).
Main search view should indicate articles that are already in the database (using checks against DOI).
Main search view should enable admin to add article to database.
From the main search view, articles are added to the table called "Main". This is the key table that contains all the artciesl from Swedish researchers.
We should probably then run a script (e.g 1 per 24 h) that retrieves all the other information for these articles and puts this information into other tables (Papers, Authors, Institutions, etc). That way, we only need 1 key table that must be correct, and from that table we populate all other tables using scripts.
When we populate the other tables (e.g. Papers) we probably need sophisticated scripts that uses chatGPT to evaluate the paper.
At a later time point: We should probably aim to download entire PubMed databas (32 million records and search for all Swedish Universities in the affiliations)


PubMex XML file format info: https://www.ncbi.nlm.nih.gov/books/NBK3828/#publisherhelp.PubMed_XML_Tagged_Format