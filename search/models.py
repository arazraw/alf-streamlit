from django.db import models

class Paper(models.Model):
    doi = models.CharField(max_length=100, primary_key=True)
    pmid = models.CharField(max_length=100)
    title = models.CharField(max_length=200)
    topic = models.CharField(max_length=100) # From Semantic Scholar
    abstract_text = models.TextField()  # Full abstract from PubMed XML
    publication_type = models.CharField(max_length=200)  # Publication type from PubMed XML
    journal_title = models.CharField(max_length=200)  # Journal title from PubMed XML
    year = models.IntegerField()  # Publication year from PubMed XML
    month = models.CharField(max_length=20)  # Publication month from PubMed XML
    citations = models.IntegerField(default=0) # From Semantic Scholar
    influential_citations = models.IntegerField(default=0) # From Semantic Scholar
    basic_science = models.BooleanField(default=False)  # Denotes a basic science paper
    trial = models.BooleanField(default=False)  # Denotes a randomized trial
    machine = models.BooleanField(default=False)  # Denotes a machine learning or AI study
    observational = models.BooleanField(default=False)  # Denotes an observational study
    guideline = models.BooleanField(default=False)  # Denotes if it is cited by guidelines

    def __str__(self):
        return self.title
    

class Author(models.Model):
    id = models.AutoField(primary_key=True) # the unique id in this specific Model
    paper = models.ForeignKey(Paper, on_delete=models.CASCADE)
    name = models.CharField(max_length=100)
    orcid_id = models.CharField(max_length=100, blank=True, null=True)  # Optional ORCID ID
    affiliation = models.CharField(max_length=200)
    citations = models.IntegerField(default=0) # From Semantic Scholar
    influential_citations = models.IntegerField(default=0) # From Semantic Scholar
    h_index = models.IntegerField(default=0)
    email = models.EmailField(blank=True, null=True)  # Optional Email

    def __str__(self):
        return self.name
