import requests
from django.core.management.base import BaseCommand
from vetu.models import AuthorImpact
from datetime import datetime
import time

class Command(BaseCommand):
    help = 'Update or create Impact instances for each paper based on Semantic Scholar data'

    BATCH_SIZE = 1000


    def fetch_author_impact_data_batch(self, s2_ids):
        api_key = 'n2iOAd5EUfimLMgOcd4B5t5wtnYLlsg9fuga7oCi'
        headers = {'X-API-KEY': api_key}

        try:
            response = requests.post(
                'https://api.semanticscholar.org/graph/v1/author/batch',
                headers=headers,
                params={'fields': 'paperCount,citationCount,hIndex'},
                json={'ids': [f"{s2_id}" for s2_id in s2_ids]}
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

    def update_author_impact_instances_batch(self, authors, data_batch):
        for s2_id, data in zip(authors, data_batch):
            if data is not None:
                try:
                    impact_instance = AuthorImpact.objects.get(s2_id=s2_id)
                    impact_instance.citations = data.get('citationCount', 0)
                    impact_instance.paper_count = data.get('paperCount', 0)
                    impact_instance.h_index = data.get('hIndex', 0)
                    impact_instance.last_updated = datetime.now()
                    impact_instance.save()
                except AuthorImpact.DoesNotExist:
                    self.stdout.write(self.style.ERROR(f'Failed to update an AuthorImpact, ensure that all author impacts exist'))
                    
            # else:
                # self.stdout.write(self.style.WARNING(f'No data returned for DOI: {paper.doi}'))

    def handle(self, *args, **options):
        authors = list(AuthorImpact.objects.all())
        total_authors = len(authors)


        completed_items = 0

        for i in range(0, total_authors, self.BATCH_SIZE):


            progress_percentage = (completed_items / total_authors) * 100
            self.stdout.write(f'\rProgress: {completed_items} out of {total_authors} complete - {progress_percentage:.2f}%', ending='')
            self.stdout.flush()


            authors_batch = authors[i:i + self.BATCH_SIZE]
            s2_ids_batch = [author.s2_id for author in authors_batch]
            data_batch = self.fetch_author_impact_data_batch(s2_ids_batch)

            self.update_author_impact_instances_batch(authors_batch, data_batch)

            # for paper, data in zip(papers_batch, data_batch):
            #     if data:
            #         self.stdout.write(self.style.SUCCESS(f'Successfully updated DOI: {paper.doi}'))
            #     else:
            #         self.stdout.write(self.style.WARNING(f'No data returned for DOI: {paper.doi}'))
            
            # Update the progress

            completed_items += self.BATCH_SIZE


            # Increment the completed items

            time.sleep(4)

        progress_percentage = 100
        self.stdout.write(f'\rProgress: {progress_percentage:.2f}% complete', ending='')
        self.stdout.flush()

        self.stdout.write(self.style.SUCCESS('\nSuccessfully updated Author Impact instances'))
