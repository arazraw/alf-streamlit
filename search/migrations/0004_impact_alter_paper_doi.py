# Generated by Django 5.0.1 on 2024-01-09 17:53

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("search", "0003_author_paper_remove_papers_main_delete_authors_and_more"),
    ]

    operations = [
        migrations.CreateModel(
            name="Impact",
            fields=[
                (
                    "paper",
                    models.OneToOneField(
                        on_delete=django.db.models.deletion.CASCADE,
                        primary_key=True,
                        serialize=False,
                        to="search.paper",
                    ),
                ),
                ("citations", models.IntegerField(default=0)),
                ("impactful_citations", models.IntegerField(default=0)),
                ("h_index", models.IntegerField(default=0)),
                ("meaningful", models.BooleanField(default=False)),
            ],
        ),
        migrations.AlterField(
            model_name="paper",
            name="doi",
            field=models.CharField(
                db_index=True, max_length=100, primary_key=True, serialize=False
            ),
        ),
    ]
