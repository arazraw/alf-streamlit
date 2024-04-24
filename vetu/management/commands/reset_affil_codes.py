from django.core.management.base import BaseCommand
from vetu.models import Paper  # Replace 'myapp' with the actual name of your Django app

class Command(BaseCommand):
    help = 'Reset affiliations_checked field for all papers'

    def handle(self, *args, **options):
        papers = Paper.objects.all()
        numb = 0

        for paper in papers:
            paper.affiliation_codes = ''
            paper.affiliations_checked = False
            paper.akademisk = False
            paper.region = False
            paper.got_un = False
            paper.kar_in = False
            paper.lun_un = False
            paper.ume_un = False
            paper.upp_un = False
            paper.lin_un = False
            paper.ore_un = False
            paper.save()
            numb += 1
            self.stdout.write(f'\rProgress: {numb}/{len(papers)} complete', ending='')
            self.stdout.flush()


        self.stdout.write(self.style.SUCCESS('\nSuccessfully reset affiliation codes field for all papers.'))
