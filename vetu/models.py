from django.db import models

class Paper(models.Model):
    doi = models.CharField(max_length=100, primary_key=True, db_index=True)
    pmid = models.CharField(max_length=100)
    title = models.CharField(max_length=200)
    topic = models.CharField(max_length=100) # From Semantic Scholar
    abstract_text = models.TextField()  # Full abstract from PubMed XML
    publication_type = models.CharField(max_length=200)  # Publication type from PubMed XML
    journal_title = models.CharField(max_length=200)  # Journal title from PubMed XML
    year = models.IntegerField()  # Publication year from PubMed XML
    month = models.CharField(max_length=20)  # Publication month from PubMed XML
    affiliations =  models.TextField()
    citations = models.IntegerField(default=0) # From Semantic Scholar
    influential_citations = models.IntegerField(default=0) # From Semantic Scholar
    basic_science = models.BooleanField(default=False)  # Denotes a basic science paper
    trial = models.BooleanField(default=False)  # Denotes a randomized trial
    machine = models.BooleanField(default=False)  # Denotes a machine learning or AI study
    observational = models.BooleanField(default=False)  # Denotes an observational study
    guideline = models.BooleanField(default=False)  # Denotes if it is cited by guidelines
    akademisk = models.BooleanField(default=False)
    region = models.BooleanField(default=False)
    got_un = models.BooleanField(default=False)
    kar_in = models.BooleanField(default=False)
    lun_un = models.BooleanField(default=False)
    ume_un = models.BooleanField(default=False)
    upp_un = models.BooleanField(default=False)
    lin_un = models.BooleanField(default=False)
    ore_un = models.BooleanField(default=False)
    affiliations_checked =  models.BooleanField(default=False)

    def __str__(self):
        return self.title
    

class Author(models.Model):
    id = models.AutoField(primary_key=True)
    paper = models.ForeignKey(Paper, on_delete=models.CASCADE, db_index=True)
    name = models.CharField(max_length=100)
    orcid_id = models.CharField(max_length=100, blank=True, null=True)  # Optional ORCID ID
    affiliation = models.CharField(max_length=200)
    email = models.EmailField(blank=True, null=True)  # Optional Email

    def __str__(self):
        return self.name


class Impact(models.Model):
    paper = models.OneToOneField(Paper, on_delete=models.CASCADE, primary_key=True)
    citations = models.IntegerField(default=0)
    impactful_citations = models.IntegerField(default=0)
    meaningful = models.BooleanField(default=False)
    citation_data = models.JSONField(default=dict, null=False)
    last_updated = models.DateTimeField(auto_now=True) 