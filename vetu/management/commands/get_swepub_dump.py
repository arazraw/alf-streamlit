from django.core.management.base import BaseCommand
from ftplib import FTP
import os
import zipfile
from django.conf import settings

class Command(BaseCommand):
    help = 'Download file from FTP server'

    def handle(self, *args, **options):
        # Display a warning and ask for confirmation
        confirm = input("THIS WILL DOWNLOAD AND EXTRACT FILES TOTALING ALMOST 30GB, ALSO THE FILE PATHS ARE HARDCODED - HAVE YOU CHECKED THEM BEFORE USE? (yes/no): ")

        if confirm.lower() != 'yes':
            self.stdout.write(self.style.ERROR("Operation aborted by user."))
            return
        else:
            self.stdout.write(self.style.NOTICE("Proceeding with download, this may take a minute."))
        

        # FTP server details
        ftp_host = 'ftp.libris.kb.se'

        # Remote file path
        remote_file_path = '/pub/spa/swepub-deduplicated.zip'

        # Local file path to save the downloaded file
        local_zip_file_path = '/run/media/Jochen/64BE-2827/swepub-deduplicated.zip'

        # Local directory path to extract the contents
        extract_dir = '/run/media/Jochen/64BE-2827/'


        ### USE THE FILEPATHS BELOW WHEN INTEGRATING INTO ACTUAL SYSTEM ON SERVER - HOWEVER DOWNLOADED FILES ARE ALL 30 GB AND EXTRACTED EVEN LARGER DONT SAVE LOCALLY
        # # Local file path to save the downloaded file
        # local_zip_file_path = os.path.join(settings.SWEPUB_DATA_DUMP_PATH, 'swepub-deduplicated.zip')

        # # Local directory path to extract the contents
        # extract_dir = os.path.join(settings.SWEPUB_DATA_DUMP_PATH, 'extracted')


        # Connect to the FTP server with anonymous access
        with FTP(ftp_host) as ftp:
            # Use anonymous login
            ftp.login()

            # Download the file
            with open(local_zip_file_path, 'wb') as local_file:
                ftp.retrbinary(f'RETR {remote_file_path}', local_file.write)

        self.stdout.write(self.style.SUCCESS(f'Successfully downloaded file from {ftp_host}:{remote_file_path} to {local_zip_file_path}.'))
    
        with zipfile.ZipFile(local_zip_file_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)

        self.stdout.write(self.style.SUCCESS(f'Successfully extracted contents to {extract_dir}.'))
