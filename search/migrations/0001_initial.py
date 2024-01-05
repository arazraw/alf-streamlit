# Generated by Django 4.2.9 on 2024-01-05 13:08

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Main',
            fields=[
                ('doi', models.CharField(max_length=100, primary_key=True, serialize=False)),
                ('pmid', models.CharField(max_length=100)),
                ('title', models.CharField(max_length=200)),
                ('year', models.IntegerField()),
            ],
        ),
        migrations.CreateModel(
            name='Papers',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('title', models.CharField(max_length=200)),
                ('topic', models.CharField(max_length=100)),
                ('citations', models.IntegerField()),
                ('main', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, to='search.main')),
            ],
        ),
        migrations.CreateModel(
            name='Authors',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=100)),
                ('affiliation', models.CharField(max_length=200)),
                ('author_citations', models.IntegerField()),
                ('main', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='search.main')),
                ('paper', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='search.papers')),
            ],
        ),
    ]