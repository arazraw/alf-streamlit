from django.db import models

class Paper(models.Model):
    doi = models.CharField(max_length=100, primary_key=True, db_index=True)
    pmid = models.CharField(max_length=100)
    title = models.CharField(max_length=500)
    topic = models.CharField(max_length=600) # Unclear
    abstract_text = models.TextField()  # Full abstract from PubMed XML
    publication_type = models.CharField(max_length=500)  # Publication type from PubMed XML
    journal_title = models.CharField(max_length=500)  # Journal title from PubMed XML
    year = models.IntegerField(null=True, blank=True, default=None)  # Publication year from PubMed XML
    month = models.CharField(max_length=20)  # Publication month from PubMed XML
    affiliations =  models.TextField()
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
    affiliation_codes = models.TextField(default='None')
    paper_fill = models.BooleanField(default=False)
    authors_checked = models.BooleanField(default=False)
    topic_codes = models.CharField(max_length=600, null=True, blank=True, default='Not Coded Yet')


    def __str__(self):
        return self.title
    

class Author(models.Model):
    id = models.AutoField(primary_key=True)
    paper = models.ForeignKey(Paper, on_delete=models.CASCADE, db_index=True)
    name = models.CharField(max_length=200)
    identifier = models.CharField(max_length=500, blank=True, null=True)  # ID saved on PubMed
    affiliation = models.CharField(max_length=10000)
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


class AuthorImpact(models.Model):
    name = models.CharField(max_length=255, null=True, blank=True)
    s2_id = models.CharField(max_length=20)
    orcid = models.CharField(max_length=20, null=True, blank=True)
    paper_count = models.IntegerField(default=0)
    citations = models.IntegerField(default=0)
    impactful_citations = models.IntegerField(default=0)
    h_index = models.IntegerField(default=0)
    last_updated = models.DateTimeField(auto_now=True) 

    def __str__(self):
        return self.s2_id #dont modify or else the update function will break