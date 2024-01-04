from django.db import models

class Main(models.Model):
    doi = models.CharField(max_length=100, primary_key=True)   # Digital Object Identifier
    pmid = models.CharField(max_length=100)  # PubMed ID
    title = models.CharField(max_length=255) # Title of the publication
    year = models.IntegerField()             # Year of publication

    def __str__(self):
        return self.title
