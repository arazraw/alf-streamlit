# update_authors.py
from django.core.management.base import BaseCommand
from vetu.models import Paper, Author
from Bio import Entrez
import requests


def fetch_authors_and_affiliations(doi):
    api_key = 'n2iOAd5EUfimLMgOcd4B5t5wtnYLlsg9fuga7oCi'
    headers = {'X-API-KEY': api_key}

    # Fetch data from Semantic Scholar
    try:
        semantic_scholar_response = requests.get(
            f'https://api.semanticscholar.org/graph/v1/paper/{doi}/authors',
            headers=headers,
            params={'fields': 'authorId,name'}
        )

        semantic_scholar_response.raise_for_status()

        semantic_scholar_data = semantic_scholar_response.json().get('data', [])
        if isinstance(semantic_scholar_data, list) and semantic_scholar_data:
            semantic_scholar_authors = semantic_scholar_data
        else:
            semantic_scholar_authors = []
            print(f'Unexpected API response for doi {doi}: {semantic_scholar_data}')
    except requests.exceptions.RequestException as err:
        semantic_scholar_authors = []
        print(f'Failed to fetch data from Semantic Scholar. Error: {err}')

    # Fetch data from PubMed
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

    pubmed_authors_data = []
    if 'AuthorList' in article:
        for author in article['AuthorList']:
            author_name = author.get('LastName', '') + ' ' + author.get('ForeName', '')
            affiliation_info = author.get('AffiliationInfo', [{}])
            affiliation = affiliation_info[0]['Affiliation'] if affiliation_info else ''
            pubmed_authors_data.append({
                'name': author_name,
                'affiliation': affiliation,
            })

    # Merge data from both sources
    merged_authors_data = []
    for semantic_scholar_author in semantic_scholar_authors:
        author_id = semantic_scholar_author.get('authorId', '')
        name = semantic_scholar_author.get('name', '')

        # Find corresponding author data from PubMed based on name
        pubmed_author_data = next((data for data in pubmed_authors_data if data['name'] == name), None)

        if pubmed_author_data:
            # Merge data from both sources
            merged_data = {
                'name': name,
                'affiliation': pubmed_author_data['affiliation'],
                'identifier': author_id
            }
            merged_authors_data.append(merged_data)
        else:
            sem_schol_data = {
                'name': name,
                'identifier': author_id,
                'affiliation': 'None'
            }
            merged_authors_data.append(sem_schol_data)

    data = {
        'ArticleTitle': title,
        'Authors': merged_authors_data
    }
    # print(data)
    return data

class Command(BaseCommand):
    help = 'Fetch authors and their affiliations for each DOI in the Paper table'

    def handle(self, *args, **options):
        total_papers = Paper.objects.filter(authors_checked=False).count()
        processed_papers = 0

        for paper in Paper.objects.filter(authors_checked=False):
            if not paper.doi:
                continue

            # Fetch authors and their affiliations for this DOI
            data = fetch_authors_and_affiliations(paper.doi)

            # Check if data is not None before entering the loop
            if data is not None:
                authors_data = data['Authors']
                self.stdout.write(f'\rProcessed {processed_papers} out of {total_papers} papers', ending='')
                self.stdout.flush()

                for author_data in authors_data:
                    # Create a new Authors object for each author
                    Author.objects.create(
                        paper=paper,
                        name=author_data['name'],
                        affiliation=author_data['affiliation'],
                        identifier=author_data['identifier']
                    )

                paper.authors_checked = True
                paper.save()
                processed_papers += 1

        self.stdout.write(f'\rProcessed {processed_papers} out of {total_papers} papers', ending='')
        self.stdout.flush()

        self.stdout.write(self.style.SUCCESS('\nSuccessfully updated authors'))
