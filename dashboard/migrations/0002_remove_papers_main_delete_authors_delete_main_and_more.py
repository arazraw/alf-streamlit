# Generated by Django 4.2.9 on 2024-01-05 18:48

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('dashboard', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='papers',
            name='main',
        ),
        migrations.DeleteModel(
            name='Authors',
        ),
        migrations.DeleteModel(
            name='Main',
        ),
        migrations.DeleteModel(
            name='Papers',
        ),
    ]