# import csv
# import re
# from django.core.management.base import BaseCommand
# from vetu.models import Paper

# class Command(BaseCommand):
#     help = 'Normalize affiliations and update affiliation_codes'

#     def handle(self, *args, **options):
#         # Load CSV file into a dictionary
#         codes_dict = {}
#         with open('/home/lukas/Documents/ALF/RegionerAkademier/affiliations_university_norm.csv', 'r') as csvfile:
#             reader = csv.reader(csvfile)
#             for row in reader:
#                 codes_dict[row[1]] = row[2]

#         # Iterate over each paper
#         numb = 0
#         paper_count = Paper.objects.count()
#         for paper in Paper.objects.all():
#             numb += 1
#             affiliations = paper.affiliations.split(';')
#             codes = set()  # Use a set to store unique codes

#             # Process each affiliation substring
#             for affiliation in affiliations:
#                 # Split affiliation on commas
#                 parts = affiliation.split(',')
#                 parts = [part.split('at') for part in parts]
#                 parts = [item for sublist in parts for item in sublist]  # Flatten the list
#                 parts = [re.sub(r'\b(?:of|the|and|och|Department|Dept|Dept.|Inst.|Institutionen|Inst)\b', '', part.lower()) for part in parts]                
#                 # Reverse iterate over parts and match codes
#                 single_digit_code_found = False
#                 digit_code = 0
#                 for part in reversed(parts):
#                     # Remove substrings containing '@' sign
#                     if '@' not in part:
#                         # Apply regex normalization
#                         normalized_part = re.sub(r'[aeiouäöü\s.\-å]', '', part.lower())
#                         # Match against codes from the CSV file
#                         for code_key, code_value in codes_dict.items():
#                             if code_key in normalized_part:
#                                 code = code_value
#                                 # Check if the first code matched is a single digit
#                                 if not single_digit_code_found and len(code) == 1:
#                                     single_digit_code_found = True
#                                     codes.add(code)  # Add code to the set
#                                     digit_code = code
#                                 # Check if all other codes in that affiliation begin with the same digit
#                                 elif single_digit_code_found and code.startswith(digit_code):
#                                     codes.add(code)  # Add code to the set
                                

#             # Update affiliation_codes field
#             paper.affiliation_codes = ';'.join(list(codes))  # Convert set to list before joining
#             paper.save()
#             self.stdout.write(f'\rProgress: {numb}/{paper_count} complete', ending='')
#             self.stdout.flush()

#         self.stdout.write(self.style.SUCCESS('Successfully updated affiliation_codes for all papers'))


import csv
import re
from django.core.management.base import BaseCommand
from vetu.models import Paper

class Command(BaseCommand):
    help = 'Normalize affiliations and update affiliation_codes'

    def handle(self, *args, **options):
        # Load CSV file into a dictionary
        codes_dict = {}
        with open('/home/lukas/Documents/ALF/RegionerAkademier/affiliations_university_norm.csv', 'r') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                codes_dict[row[1]] = row[2]

        # Iterate over each paper
        numb = 0
        paper_count = Paper.objects.count()
        for paper in Paper.objects.all():
            numb += 1
            affiliations = paper.affiliations.split(';')
            codes = set()  # Use a set to store unique codes
            # Process each affiliation substring
            for affiliation in affiliations:
                # Split affiliation on commas
                parts = affiliation.split(',')
                parts = [part.split('at') for part in parts]
                parts = [item for sublist in parts for item in sublist]  # Flatten the list
                parts = [re.sub(r'\b(?:of|the|and|och|Department|Dept|Dept.|Inst.|Institutionen|Inst)\b', '', part.lower()) for part in parts]                
                # Reverse iterate over parts and match codes
                single_digit_code_found = False
                digit_code = 0
                temp_codes = set()
                for part in reversed(parts):
                    # Remove substrings containing '@' sign
                    if '@' not in part:
                        # Apply regex normalization
                        normalized_part = re.sub(r'[aeiouäöü\s.\-å]', '', part.lower())
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
            paper.affiliation_codes = ';'.join(list(codes))  # Convert set to list before joining
            paper.save()
            self.stdout.write(f'\rProgress: {numb}/{paper_count} complete', ending='')
            self.stdout.flush()

        self.stdout.write(self.style.SUCCESS('Successfully updated affiliation_codes for all papers'))
