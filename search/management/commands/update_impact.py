import requests
from django.core.management.base import BaseCommand
from search.models import Paper
import time


class Command(BaseCommand):
    help = 'Update citation count, influential citation count, and topic for each paper'

    def handle(self, *args, **options):
        api_key = 'n2iOAd5EUfimLMgOcd4B5t5wtnYLlsg9fuga7oCi'
        headers = {'X-API-KEY': api_key}

        for paper in Paper.objects.all():
            doi = paper.doi
            if doi:
                try:
                    response = requests.post(
                        'https://api.semanticscholar.org/graph/v1/paper/batch',
                        headers=headers,
                        params={'fields': 'citationCount,influentialCitationCount'},  # Removed 'topics' from fields
                        json={'ids': [f"DOI:{doi}"]}
                    )

                    if response.status_code != 200:
                        self.stdout.write(self.style.ERROR(f'Response Content for DOI {doi}: {response.text}'))
                        continue

                    response_data = response.json()

                    if isinstance(response_data, list) and response_data:
                        data = response_data[0]
                        paper.citations = data.get('citationCount', 0)
                        paper.influential_citations = data.get('influentialCitationCount', 0)
                        paper.save()
                    else:
                        self.stdout.write(self.style.WARNING(f'No data returned for DOI: {doi}'))
                except requests.exceptions.RequestException as err:
                    self.stdout.write(self.style.ERROR(f'Failed to fetch data for DOI: {doi}. Error: {err}'))
            else:
                self.stdout.write(self.style.WARNING(f'Skipping paper with DOI: None'))
            
            # Sleep for 3 seconds
            time.sleep(3)
