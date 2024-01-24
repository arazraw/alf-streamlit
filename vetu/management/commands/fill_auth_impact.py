from django.core.management.base import BaseCommand
from vetu.models import Author, AuthorImpact

class Command(BaseCommand):
    help = 'Create AuthorImpact instances for each unique identifier in Author'

    def handle(self, *args, **options):
        # Get all unique identifiers and names from the Author model
        authors_data = Author.objects.values('identifier', 'name').distinct()

        for author_data in authors_data:
            identifier = author_data['identifier']
            name = author_data['name']

            # Check if an AuthorImpact instance with the same identifier already exists
            if not AuthorImpact.objects.filter(s2_id=identifier).exists():
                # Ensure identifier is not None before creating an AuthorImpact instance
                if identifier is not None:
                    # If not, create a new AuthorImpact instance
                    author_impact = AuthorImpact(s2_id=identifier, name=name)
                    author_impact.full_clean()  # Validate the model
                    author_impact.save()

                    self.stdout.write(self.style.SUCCESS(f'Successfully created AuthorImpact for identifier {identifier}'))
                else:
                    self.stdout.write(self.style.WARNING(f'Skipping AuthorImpact creation for None identifier'))
            else:
                self.stdout.write(self.style.WARNING(f'AuthorImpact for identifier {identifier} already exists'))
