import requests
from django.core.management.base import BaseCommand
from search.models import Paper, Impact
from datetime import datetime
import time

class Command(BaseCommand):
    help = 'Update or create Impact instances for each paper based on Semantic Scholar data'

    def fetch_paper_data(self, doi):
        api_key = 'n2iOAd5EUfimLMgOcd4B5t5wtnYLlsg9fuga7oCi'
        headers = {'X-API-KEY': api_key}

        try:
            response = requests.post(
                'https://api.semanticscholar.org/graph/v1/paper/batch',
                headers=headers,
                params={'fields': 'citationCount,influentialCitationCount,citations'},  # Adjust the fields as needed
                json={'ids': [f"DOI:{doi}"]}
            )

            response.raise_for_status()

            response_data = response.json()
            if isinstance(response_data, list) and response_data:
                return response_data[0]
            else:
                self.stdout.write(self.style.WARNING(f'Unexpected API response for DOI {doi}: {response_data}'))
        except requests.exceptions.RequestException as err:
            self.stdout.write(self.style.ERROR(f'Failed to fetch data for DOI: {doi}. Error: {err}'))

        return None

    def create_or_update_impact_instance(self, paper, data):
        try:
            impact_instance = Impact.objects.get(paper=paper)
            impact_instance.citations = data.get('citationCount', 0)
            impact_instance.impactful_citations = data.get('influentialCitationCount', 0)
            impact_instance.meaningful = False
            impact_instance.citation_data = data.get('citations', [])
            impact_instance.last_updated = datetime.now()  # Add the last_updated field
            impact_instance.save()
        except Impact.DoesNotExist:
            Impact.objects.create(
                paper=paper,
                citations=data.get('citationCount', 0),
                impactful_citations=data.get('influentialCitationCount', 0),
                meaningful=False,
                citation_data=data.get('citations', []),
                last_updated=datetime.now()  # Add the last_updated field
            )

    def handle(self, *args, **options):
        for paper in Paper.objects.exclude(doi__isnull=True).exclude(doi__exact='').exclude(doi__iexact='None'):
            doi = paper.doi
            data = self.fetch_paper_data(doi)
            if data is not None:
                self.create_or_update_impact_instance(paper, data)
                self.stdout.write(self.style.SUCCESS(f'Successfully updated DOI: {doi}'))
            else:
                self.stdout.write(self.style.WARNING(f'No data returned for DOI: {doi}'))

            time.sleep(4)

        self.stdout.write(self.style.SUCCESS('Successfully created or updated impact instances'))
