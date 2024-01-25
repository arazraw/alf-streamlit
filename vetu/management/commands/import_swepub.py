# import_swepub.py

from django.core.management.base import BaseCommand
import json
from pathlib import Path
from django.conf import settings
from vetu.models import Paper
from dateutil import parser  # Import the dateutil library for parsing dates

class Command(BaseCommand):
    help = 'Extract information from each entry of the .jsonl file'

    def extract_year_from_date(self, date_str):
        try:
            parsed_date = parser.parse(date_str)
            return parsed_date.year
        except ValueError:
            # Handle the case where parsing fails, return None or handle accordingly
            return None

    def extract_topic_codes(self, subject_data_list):
        # Extracts all topic codes from the subject data list and concatenates them into a string
        topic_codes = []
        for subject in subject_data_list:
            if 'Topic' in subject.get('@type', ''):
                topic_code = subject.get('code', '')
                if topic_code:
                    topic_codes.append(topic_code)
        return ','.join(topic_codes)

    def handle(self, *args, **options):
        # Local directory path to the extracted .jsonl file
        jsonl_file_path = Path('/run/media/Jochen/Elements/swepub/') / 'swepub-deduplicated.jsonl'

        try:
            with open(jsonl_file_path, 'r', encoding='utf-8') as jsonl_file:
                # Loop through each entry in the .jsonl file
                for line in jsonl_file:
                    entry = json.loads(line)

                    # Extract information from the nested structure
                    master_data = entry.get('master', {})
                    instanceOf = master_data.get('instanceOf', {})
                    subject_data_list = instanceOf.get('subject', [])

                    # Extract topic codes as a string
                    topic_codes = self.extract_topic_codes(subject_data_list)

                    # Extract specific fields
                    doi = ''
                    pmid = ''
                    year = ''

                    identified_by_list = master_data.get('identifiedBy', [])

                    for identifier in identified_by_list:
                        if identifier.get('@type', '') == 'DOI':
                            doi = identifier.get('value', '')
                            doi = doi.replace('https://doi.org/', '')
                        elif identifier.get('@type', '') == 'PMID':
                            pmid = identifier.get('value', '')

                    publication_list = master_data.get('publication', [])

                    for identifier in publication_list:
                        if identifier.get('@type', '') == 'Publication':
                            year = identifier.get('date', '')

                    # Extract the year from the date
                    year = self.extract_year_from_date(year)

                    # Check if either DOI or PMID exists
                    if not doi and not pmid:
                        self.stdout.write(self.style.SUCCESS('Both DOI and PMID are non-existent. Skipping entry.'))
                        continue

                    # Display the extracted information
                    self.stdout.write(self.style.SUCCESS(f'DOI: {doi}'))
                    self.stdout.write(self.style.SUCCESS(f'PMID: {pmid}'))
                    self.stdout.write(self.style.SUCCESS(f'Year: {year}'))
                    self.stdout.write(self.style.SUCCESS(f'Topic Codes: {topic_codes}'))

                    paper_exists = Paper.objects.filter(doi=doi) or Paper.objects.filter(pmid=pmid)

                    if not paper_exists:
                        # Create a new entry in the Paper model
                        new_paper = Paper(doi=doi, pmid=pmid, year=year, topic_codes=topic_codes)
                        new_paper.save()
                        self.stdout.write(self.style.SUCCESS(f'New entry created in the Paper model.'))
                    else:
                        self.stdout.write(self.style.SUCCESS('DOI or PMID already exists in the Paper model.'))

        except FileNotFoundError:
            self.stdout.write(self.style.ERROR(f'Error: File not found - {jsonl_file_path}'))

        except json.JSONDecodeError:
            self.stdout.write(self.style.ERROR(f'Error: Invalid JSON format in {jsonl_file_path}'))

        except Exception as e:
            self.stdout.write(self.style.ERROR(f'Error: {str(e)}'))
