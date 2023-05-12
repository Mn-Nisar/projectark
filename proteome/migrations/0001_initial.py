# Generated by Django 4.1.4 on 2023-05-12 06:45

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Contaminant',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(default='cont', max_length=255)),
                ('file', models.FileField(null=True, upload_to='contaminant/')),
            ],
        ),
        migrations.CreateModel(
            name='DataAnalysis',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('file', models.FileField(null=True, upload_to='documents/')),
                ('created', models.DateTimeField(auto_now=True)),
                ('user', models.CharField(max_length=255, null=True)),
                ('resultData', models.FileField(null=True, upload_to='documents/')),
                ('expression_data', models.FileField(null=True, upload_to='documents/')),
                ('labledData', models.BooleanField(default=False)),
                ('lableFreeData', models.BooleanField(default=False)),
                ('geneOntology', models.FileField(null=True, upload_to='documents/')),
                ('stringdb', models.FileField(null=True, upload_to='documents/')),
                ('clusterdb', models.FileField(null=True, upload_to='documents/')),
                ('ibaq', models.FileField(null=True, upload_to='documents/')),
                ('fastafile', models.FileField(null=True, upload_to='documents/')),
                ('dleteted_resd', models.FileField(null=True, upload_to='documents/')),
            ],
        ),
        migrations.CreateModel(
            name='Example',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255, null=True)),
                ('file', models.FileField(null=True, upload_to='example/')),
                ('usethis', models.BooleanField(default=True)),
                ('labledData', models.BooleanField(default=True)),
                ('lableFreeData', models.BooleanField(default=False)),
                ('number_of_sample', models.IntegerField(default=0)),
                ('biological_rep', models.BooleanField(default=True)),
                ('number_of_batches', models.IntegerField(default=0)),
                ('fastafile', models.FileField(blank=True, null=True, upload_to='example/')),
            ],
        ),
        migrations.CreateModel(
            name='FastaBackup',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('species', models.CharField(max_length=100, null=True)),
                ('fasta_file', models.FileField(null=True, upload_to='fasta/')),
            ],
        ),
        migrations.CreateModel(
            name='Forzip',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('analysis_id', models.BigIntegerField()),
                ('files', models.FileField(upload_to='documents/')),
            ],
        ),
        migrations.CreateModel(
            name='PlotExample',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('plot_type', models.CharField(max_length=255, null=True)),
                ('file', models.FileField(null=True, upload_to='example/')),
            ],
        ),
        migrations.CreateModel(
            name='Ploting',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('plot_type', models.CharField(max_length=100, null=True)),
                ('file', models.FileField(null=True, upload_to='documents/')),
                ('user', models.CharField(max_length=255, null=True)),
            ],
        ),
    ]