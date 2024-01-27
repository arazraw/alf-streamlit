from django.core.management.base import BaseCommand
from vetu.models import Paper

class Command(BaseCommand):
    help = 'Handle topics in the Paper model'

    def handle(self, *args, **options):
        # Retrieve all papers
        papers = Paper.objects.all()

        for paper in papers:
            # Retrieve topics for the current paper
            topics = paper.topic_codes

            # Split topics by ',' and remove duplicates
            if topics:
                unique_topics = sorted(set(topics.split(',')))

                # Join unique topics with ';' separator
                new_topics = ';'.join(unique_topics)

                # Update topics for the current paper
                paper.topic_codes = new_topics
                paper.save()

        self.stdout.write(self.style.SUCCESS('Topics handling completed'))
