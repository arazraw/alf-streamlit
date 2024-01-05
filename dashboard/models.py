from django.db import models

class Main(models.Model):
    doi = models.CharField(max_length=100, primary_key=True)
    pmid = models.CharField(max_length=100)
    title = models.CharField(max_length=200)
    year = models.IntegerField()

    def __str__(self):
        return self.title

class Papers(models.Model):
    id = models.AutoField(primary_key=True)
    main = models.OneToOneField(Main, on_delete=models.CASCADE)
    title = models.CharField(max_length=200)
    topic = models.CharField(max_length=100)
    citations = models.IntegerField()

    def __str__(self):
        return self.title

class Authors(models.Model):
    id = models.AutoField(primary_key=True)
    main = models.ForeignKey(Main, on_delete=models.CASCADE)
    paper = models.ForeignKey(Papers, on_delete=models.CASCADE)
    name = models.CharField(max_length=100)
    affiliation = models.CharField(max_length=200)
    author_citations = models.IntegerField()

    def __str__(self):
        return self.name
