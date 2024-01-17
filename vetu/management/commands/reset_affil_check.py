from django.core.management.base import BaseCommand
from vetu.models import Paper  # Replace 'myapp' with the actual name of your Django app

class Command(BaseCommand):
    help = 'Reset affiliations_checked field for all papers'

    def handle(self, *args, **options):
        papers = Paper.objects.all()

        for paper in papers:
            paper.affiliations_checked = False
            paper.save()

        self.stdout.write(self.style.SUCCESS('Successfully reset affiliations_checked field for all papers.'))
