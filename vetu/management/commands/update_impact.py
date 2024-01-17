import requests
from django.core.management.base import BaseCommand
from vetu.models import Paper, Impact
from datetime import datetime
import time

class Command(BaseCommand):
    help = 'Update or create Impact instances for each paper based on Semantic Scholar data'

    BATCH_SIZE = 100

    def fetch_paper_data_batch(self, dois):
        api_key = 'n2iOAd5EUfimLMgOcd4B5t5wtnYLlsg9fuga7oCi'
        headers = {'X-API-KEY': api_key}

        try:
            response = requests.post(
                'https://api.semanticscholar.org/graph/v1/paper/batch',
                headers=headers,
                params={'fields': 'citationCount,influentialCitationCount,citations'},
                json={'ids': [f"DOI:{doi}" for doi in dois]}
            )

            response.raise_for_status()

            response_data = response.json()
            if isinstance(response_data, list) and response_data:
                return response_data
            else:
                self.stdout.write(self.style.WARNING(f'Unexpected API response for batch: {response_data}'))
        except requests.exceptions.RequestException as err:
            self.stdout.write(self.style.ERROR(f'Failed to fetch data for batch. Error: {err}'))

        return []

    def create_or_update_impact_instances_batch(self, papers, data_batch):
        for paper, data in zip(papers, data_batch):
            try:
                impact_instance = Impact.objects.get(paper=paper)
                impact_instance.citations = data.get('citationCount', 0)
                impact_instance.impactful_citations = data.get('influentialCitationCount', 0)
                impact_instance.meaningful = False
                impact_instance.citation_data = data.get('citations', [])
                impact_instance.last_updated = datetime.now()
                impact_instance.save()
            except Impact.DoesNotExist:
                Impact.objects.create(
                    paper=paper,
                    citations=data.get('citationCount', 0),
                    impactful_citations=data.get('influentialCitationCount', 0),
                    meaningful=False,
                    citation_data=data.get('citations', []),
                    last_updated=datetime.now()
                )

    def handle(self, *args, **options):
        papers = list(Paper.objects.exclude(doi__isnull=True).exclude(doi__exact='').exclude(doi__iexact='None'))
        total_papers = len(papers)

        for i in range(0, total_papers, self.BATCH_SIZE):
            papers_batch = papers[i:i + self.BATCH_SIZE]
            dois_batch = [paper.doi for paper in papers_batch]
            data_batch = self.fetch_paper_data_batch(dois_batch)

            self.create_or_update_impact_instances_batch(papers_batch, data_batch)

            for paper, data in zip(papers_batch, data_batch):
                if data:
                    self.stdout.write(self.style.SUCCESS(f'Successfully updated DOI: {paper.doi}'))
                else:
                    self.stdout.write(self.style.WARNING(f'No data returned for DOI: {paper.doi}'))

            time.sleep(4)

        self.stdout.write(self.style.SUCCESS('Successfully created or updated impact instances'))
