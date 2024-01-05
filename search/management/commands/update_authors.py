# update_authors.py
from django.core.management.base import BaseCommand
from search.models import Main, Authors
from Bio import Entrez


def fetch_authors_and_affiliations(doi):
    handle = Entrez.esearch(db="pubmed", term=doi)
    record = Entrez.read(handle)
    handle.close()

    id_list = record.get("IdList", [])
    if not id_list:
        print(f"No articles found for DOI '{doi}'.")
        return None

    handle = Entrez.efetch(db="pubmed", id=id_list[0], rettype="medline", retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    record = records['PubmedArticle'][0]
    article = record['MedlineCitation']['Article']
    title = article.get('ArticleTitle', 'No title available')

    authors_data = []
    if 'AuthorList' in article:
        for author in article['AuthorList']:
            author_name = author.get('LastName', '') + ' ' + author.get('ForeName', '')
            affiliation = author.get('AffiliationInfo', [{}])[0].get('Affiliation', '')
            authors_data.append({
                'name': author_name,
                'affiliation': affiliation,
                'author_citations': 0  # Replace with actual citation count if available
            })

    data = {
        'ArticleTitle': title,
        'Authors': authors_data,
    }

    return data

class Command(BaseCommand):
    help = 'Fetch authors and their affiliations for each DOI in the Main table'

    def handle(self, *args, **options):
        for main in Main.objects.all():
            # Fetch authors and their affiliations for this DOI
            authors_data = fetch_authors_and_affiliations(main.doi)

            for author_data in authors_data['Authors']:
                # Create a new Authors object for each author
                Authors.objects.create(
                    main=main,
                    paper=main.papers,  # Replace with actual paper
                    name=author_data['name'],
                    affiliation=author_data['affiliation'],
                    author_citations=author_data['author_citations']
                )

        self.stdout.write(self.style.SUCCESS('Successfully updated authors'))