from django.core.management.base import BaseCommand
from search.models import Main, Papers  # Update this line
from Bio import Entrez


def fetch_article_by_doi(doi):
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

    data = {
        'ArticleTitle': title,
    }

    return data



class Command(BaseCommand):
    help = 'Updates the Papers table with missing DOIs from the Main table'

    def handle(self, *args, **options):
        for main in Main.objects.all():
            if not Papers.objects.filter(main=main).exists():
                data = fetch_article_by_doi(main.doi)
                Papers.objects.create(main=main, title=data['ArticleTitle'])

        self.stdout.write(self.style.SUCCESS('Successfully updated Papers table'))