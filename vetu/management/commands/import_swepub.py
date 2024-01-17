# import_swepub.py

from django.core.management.base import BaseCommand
import json
from pathlib import Path
from django.conf import settings

class Command(BaseCommand):
    help = 'Load and display the first entry from the extracted .jsonl file'

    def handle(self, *args, **options):
        # Local directory path to the extracted .jsonl file
        jsonl_file_path = Path('/run/media/Jochen/Elements/swepub/')/ 'swepub-deduplicated.jsonl'
        # jsonl_file_path = Path(settings.SWEPUB_DATA_DUMP_PATH) / 'extracted' / 'swepub-deduplicated.jsonl'
        output_txt_file_path = '/run/media/Jochen/Elements/swepub/output.txt'

        try:
            with open(jsonl_file_path, 'r', encoding='utf-8') as jsonl_file:
                # Read the first line from the .jsonl file
                first_entry = json.loads(jsonl_file.readline())

                # # Display the first entry
                # self.stdout.write(self.style.SUCCESS('First entry from .jsonl file:'))
                # self.stdout.write(self.style.SUCCESS(json.dumps(first_entry, indent=2)))

                with open(output_txt_file_path, 'w', encoding='utf-8') as output_txt_file:
                    output_txt_file.write(json.dumps(first_entry, indent=2))

                self.stdout.write(self.style.SUCCESS(f'First entry written to {output_txt_file_path}.'))


        except FileNotFoundError:
            self.stdout.write(self.style.ERROR(f'Error: File not found - {jsonl_file_path}'))

        except json.JSONDecodeError:
            self.stdout.write(self.style.ERROR(f'Error: Invalid JSON format in {jsonl_file_path}'))

        except Exception as e:
            self.stdout.write(self.style.ERROR(f'Error: {str(e)}'))
