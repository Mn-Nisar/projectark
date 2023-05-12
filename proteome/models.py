from django.db import models
class DataAnalysis(models.Model):
    file = models.FileField(null = True, upload_to='documents/')
    created = models.DateTimeField(auto_now=True)
    user = models.CharField(null = True, max_length = 255)
    resultData = models.FileField(null = True, upload_to='documents/')
    expression_data = models.FileField(null = True, upload_to='documents/')
    labledData = models.BooleanField(default=False)
    lableFreeData = models.BooleanField(default=False)
    geneOntology =  models.FileField(null = True, upload_to='documents/')
    stringdb = models.FileField(null = True, upload_to = 'documents/')
    clusterdb = models.FileField(null = True, upload_to = 'documents/')
    ibaq = models.FileField(null = True, upload_to = 'documents/')
    fastafile = models.FileField(null = True, upload_to = 'documents/')
    dleteted_resd = models.FileField(null = True, upload_to = 'documents/')

    def __str__(self):
        return str(self.user)

class Example(models.Model):
    name = models.CharField(null = True, max_length = 255)
    file = models.FileField(null = True, upload_to='example/')
    usethis = models.BooleanField(default=True)
    labledData = models.BooleanField(default=True)
    lableFreeData = models.BooleanField(default=False)
    number_of_sample = models.IntegerField(default = 0)
    biological_rep = models.BooleanField(default=True)
    number_of_batches = models.IntegerField(default = 0)
    fastafile = models.FileField(null = True, blank = True, upload_to='example/')

    def __str__(self):
        if self.name:
            return self.name
        else:
            return None


class Contaminant(models.Model):
    name = models.CharField(default = "cont", max_length = 255)
    file = models.FileField(null = True, upload_to='contaminant/')


    def __str__(self):
        if self.name:
            return self.name
        else:
            return None

class FastaBackup(models.Model):
    species = models.CharField(max_length = 100, null = True)
    fasta_file = models.FileField(null = True, upload_to='fasta/')

    def __str__(self):
        if self.species:
            return self.species
        else:
            return None


class Ploting(models.Model):
    plot_type = models.CharField(max_length= 100, null=True)
    file = models.FileField(null = True, upload_to='documents/')
    user = models.CharField(null = True, max_length = 255)

    def __str__(self):
        if self.plot_type:
            return self.plot_type
        else:
            return None


class Forzip(models.Model):
    analysis_id = models.BigIntegerField()
    files = models.FileField(upload_to='documents/')

    def __str__(self):
        if self.analysis_id:
            return str(self.analysis_id)
        else:
            return None


class PlotExample(models.Model):
    plot_type = models.CharField(null = True, max_length = 255)
    file = models.FileField(null = True, upload_to='example/')

    def __str__(self):
        if self.plot_type:
            return self.plot_type
        else:
            return 'None'
        