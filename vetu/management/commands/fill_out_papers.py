# fill out papers.py
from django.core.management.base import BaseCommand
from vetu.models import Paper
from Bio import Entrez
from urllib.error import HTTPError

def fetch_paper_information(identifier, identifier_type):
    try:
        if identifier_type == 'DOI':
            handle = Entrez.esearch(db="pubmed", term=identifier)        
            record = Entrez.read(handle)
            handle.close()

            id_list = record.get("IdList", [])
            if not id_list:
                print(f"No articles found for {identifier_type} '{identifier}'.")
                return None

            handle = Entrez.efetch(db="pubmed", id=id_list[0], rettype="medline", retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            record = records['PubmedArticle'][0]

            doi = identifier
            pmid = record['MedlineCitation']['PMID']
            
        elif identifier_type == 'PMID':
            handle = Entrez.esearch(db="pubmed", term=identifier, field="pmid")
            record = Entrez.read(handle)
            handle.close()

            id_list = record.get("IdList", [])
            if not id_list:
                print(f"No articles found for {identifier_type} '{identifier}'.")
                return None

            handle = Entrez.efetch(db="pubmed", id=id_list[0], rettype="medline", retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            record = records['PubmedArticle'][0]
            article = record['MedlineCitation']['Article']

            pmid = identifier
            doi = None
            if 'ELocationID' in article:
                for location in article['ELocationID']:
                    if isinstance(location, str) and location.startswith('10.'):
                        doi = location
                        break

        else:
            return None


        
        article = record['MedlineCitation']['Article']
        title = article.get('ArticleTitle', 'No title available')
        abstract = article.get('Abstract', {}).get('AbstractText', [''])[0]
        publication_type = article.get('PublicationTypeList', [''])[0]
        journal_title = article.get('Journal', {}).get('Title', 'No journal title available')
        publication_date = article.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
        month = publication_date.get('Month', '')
        # Extracting all affiliations for all authors
        all_affiliations = []
        if 'AuthorList' in article:
            for author in article['AuthorList']:
                if 'AffiliationInfo' in author:
                    author_affiliations = [info['Affiliation'] for info in author['AffiliationInfo']]
                    all_affiliations.extend(author_affiliations)

        unique_affiliations = list(set(all_affiliations))

        data = {
            'ArticleTitle': title,
            'DOI': doi,
            'PMID': pmid,
            'Abstract': abstract,
            'PublicationType': publication_type,
            'JournalTitle': journal_title,
            'Month': month,
            'Affiliations':  '; '.join(unique_affiliations)
            # Add more fields as needed
        }

        

        return data
    
    except HTTPError as e:
        print(f"HTTP Error: {e.code} - {e.reason}. Unable to fetch information for {identifier_type} '{identifier}'.")
        return None

class Command(BaseCommand):
    help = 'Fetch additional information for each DOI in the Paper table'

    def handle(self, *args, **options):
        for paper in Paper.objects.filter(title=''):

            self.stdout.write(self.style.NOTICE(f'Requesting PubMed for Id {paper.doi} {paper.pmid}'))
            # Fetch additional information for this DOI
            if paper.doi:
                data = fetch_paper_information(paper.doi, 'DOI')
            # If DOI is missing, try searching with PMID
            elif paper.pmid:
                data = fetch_paper_information(paper.pmid, 'PMID')


            # Check if data is not None before updating the Paper object
            if data is not None:
                paper.title = data['ArticleTitle']
                paper.abstract_text = data['Abstract']
                paper.publication_type = data['PublicationType']
                paper.journal_title = data['JournalTitle']
                paper.month = data['Month']
                paper.affiliations = data['Affiliations']
                paper.doi = data['DOI']
                paper.pmid = data['PMID']
                # Update more fields as needed
                paper.save()


        self.stdout.write(self.style.SUCCESS('Successfully updated Paper objects'))
