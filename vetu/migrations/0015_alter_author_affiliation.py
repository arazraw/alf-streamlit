# Generated by Django 4.2.9 on 2024-02-05 10:34

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("vetu", "0014_alter_author_affiliation"),
    ]

    operations = [
        migrations.AlterField(
            model_name="author",
            name="affiliation",
            field=models.CharField(max_length=10000),
        ),
    ]
