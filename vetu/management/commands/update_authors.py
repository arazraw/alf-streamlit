# update_authors.py
from django.core.management.base import BaseCommand
from vetu.models import Paper, Author
from Bio import Entrez

def fetch_authors_and_affiliations(doi):
    print(f"Fetching data for DOI: {doi}")
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
            affiliation_info = author.get('AffiliationInfo', [{}])
            affiliation = affiliation_info[0]['Affiliation'] if affiliation_info else ''
            authors_data.append({
                'name': author_name,
                'affiliation': affiliation
            })

    data = {
        'ArticleTitle': title,
        'Authors': authors_data,
    }

    return data

class Command(BaseCommand):
    help = 'Fetch authors and their affiliations for each DOI in the Paper table'

    def handle(self, *args, **options):
        for paper in Paper.objects.all():
            if not paper.doi:
                print(f"Invalid DOI for paper {paper.title}")
                continue
            # Fetch authors and their affiliations for this DOI
            data = fetch_authors_and_affiliations(paper.doi)

            # Check if data is not None before entering the loop
            if data is not None:
                authors_data = data['Authors']
                for author_data in authors_data:
                    # Create a new Authors object for each author
                    Author.objects.create(
                        paper=paper,
                        name=author_data['name'],
                        affiliation=author_data['affiliation']
                    )

        self.stdout.write(self.style.SUCCESS('Successfully updated authors'))