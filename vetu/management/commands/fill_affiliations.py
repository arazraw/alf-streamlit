from django.core.management.base import BaseCommand
from vetu.models import Paper 

class Command(BaseCommand):
    help = 'Update affiliation flags based on affiliations field'

    def handle(self, *args, **options):
        papers = Paper.objects.filter(affiliations_checked=False)

        for paper in papers:
            affiliations = paper.affiliations.lower()

            # Check for specific strings in the affiliations field
            # When adding fields to this and wishing to rerun over all in Paper there is a command called reset_affil_check
            paper.akademisk = any(keyword in affiliations for keyword in ['akademisk', 'akademi', 'university'])
            paper.region = any(keyword in affiliations for keyword in ['region', 'hospital'])
            paper.got_un = any(keyword in affiliations for keyword in ['university of gothenburg', 'gothenburg university', 'göteborgs universitet'])
            paper.kar_in = any(keyword in affiliations for keyword in ['karolinska institutet', 'karolinska institute'])
            paper.lun_un = any(keyword in affiliations for keyword in ['lund universitet', 'lund university'])
            paper.ume_un = any(keyword in affiliations for keyword in ['umeå university', 'umea university'])
            paper.upp_un = any(keyword in affiliations for keyword in ['uppsala university', 'uppsala universitet'])
            paper.lin_un = any(keyword in affiliations for keyword in ['linköping university', 'linköping universitet'])
            paper.ore_un = any(keyword in affiliations for keyword in ['örebro university', 'örebro universitet'])

            # Set affiliations_checked to True to indicate that the affiliations field has been processed
            paper.affiliations_checked = True

            # Save the changes
            paper.save()

        self.stdout.write(self.style.SUCCESS('Successfully updated affiliation flags.'))