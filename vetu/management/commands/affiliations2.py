import csv
import re
from django.core.management.base import BaseCommand
from vetu.models import Paper
from django.conf import settings
import os

class Command(BaseCommand):
    help = 'Normalize affiliations and update affiliation_codes'

    def handle(self, *args, **options):
        # Load CSV file into a dictionary
        codes_dict = {}
        file_path = os.path.join(settings.BASE_DIR, 'RegionerAkademier/affiliations_university_norm.csv')

        with open(file_path, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                codes_dict[row[1]] = row[2]

        # Iterate over each paper
        numb = 0
        paper_count = Paper.objects.count()
        papers = Paper.objects.filter(affiliations_checked=False)
        for paper in papers:
            numb += 1
            affiliations = paper.affiliations.split(';')
            codes = set()  # Use a set to store unique codes
            # Process each affiliation substring
            for affiliation in affiliations:
                # Split affiliation on commas
                parts = affiliation.split(',')
                parts = [part.split('at') for part in parts]
                parts = [item for sublist in parts for item in sublist]  # Flatten the list
                parts = [re.sub(r'\b(?:of|the|and|och|department|dept|dept\.|inst\.|institutionen|inst|för)\b', '', part.lower()) for part in parts]                
                # Reverse iterate over parts and match codes
                single_digit_code_found = False
                digit_code = 0
                temp_codes = set()
                for part in reversed(parts):
                    # Remove substrings containing '@' sign
                    if '@' not in part:
                        # Apply regex normalization
                        input_string = re.sub(r'\(.*?\)', '', part)  # Remove parentheses and anything inside
                        input_string = re.sub(r'[-,]', '', input_string)  # Remove hyphens and commas
                        normalized_part = re.sub(r'[aeiouäöü\s.\-"\'å]', '', input_string.lower())
                        # Match against codes from the CSV file
                        for code_key, code_value in codes_dict.items():
                            if code_key in normalized_part:
                                code = code_value
                                # Check if the code is a single digit
                                if not single_digit_code_found:
                                    if len(code) == 1:
                                        single_digit_code_found = True
                                        codes.add(code)  # Add single digit code to the set
                                        digit_code = code
                                        # Add all codes in temp_codes that match the first digit
                                        for temp_code in temp_codes:
                                            if temp_code.startswith(digit_code):
                                                codes.add(temp_code)
                                        # Reset temp_codes
                                        temp_codes = set()
                                    else:
                                        temp_codes.add(code)  # Add code to temp set                                
                                # Check if the code matches the first digit of the single digit code found
                                elif single_digit_code_found and code.startswith(digit_code):
                                    codes.add(code)  # Add code to the set
                                

            # Update affiliation_codes field
            affiliations_c = paper.affiliations.lower()

            paper.affiliation_codes = ';'.join(list(codes))  # Convert set to list before joining
            paper.akademisk = any(keyword in affiliations_c for keyword in ['akademisk', 'akademi', 'university', 'högskola'])
            paper.region = any(keyword in affiliations_c for keyword in ['region', 'hospital'])
            paper.got_un = any(keyword in affiliations_c for keyword in ['university of gothenburg', 'gothenburg university', 'göteborgs universitet'])
            paper.kar_in = any(keyword in affiliations_c for keyword in ['karolinska institutet', 'karolinska institute'])
            paper.lun_un = any(keyword in affiliations_c for keyword in ['lund universitet', 'lund university'])
            paper.ume_un = any(keyword in affiliations_c for keyword in ['umeå university', 'umea university'])
            paper.upp_un = any(keyword in affiliations_c for keyword in ['uppsala university', 'uppsala universitet'])
            paper.lin_un = any(keyword in affiliations_c for keyword in ['linköping university', 'linköping universitet'])
            paper.ore_un = any(keyword in affiliations_c for keyword in ['örebro university', 'örebro universitet'])
            single_digit_codes = [code for code in codes if len(code) == 1]
            for code in single_digit_codes:
                if code == '1':
                    paper.got_un = True
                    paper.akademisk = True
                elif code == '2':
                    paper.kar_in = True
                    paper.akademisk = True
                elif code == '3':
                    paper.lin_un = True
                    paper.akademisk = True
                elif code == '4':
                    paper.ume_un = True
                    paper.akademisk = True
                elif code == '5':
                    paper.ore_un = True
                    paper.akademisk = True
                elif code == '6':
                    paper.upp_un = True
                    paper.akademisk = True
                elif code == '7':
                    paper.lun_un = True
                    paper.akademisk = True
            # Set affiliations_checked to True to indicate that the affiliations field has been processed
            paper.affiliations_checked = True

            paper.save()
            self.stdout.write(f'\rProgress: {numb}/{paper_count} complete', ending='')
            self.stdout.flush()

        self.stdout.write(self.style.SUCCESS('\nSuccessfully updated affiliation_codes for all papers'))
