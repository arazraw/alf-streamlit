# Generated by Django 5.0.1 on 2024-01-27 20:28

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ("vetu", "0012_paper_topic_codes"),
    ]

    operations = [
        migrations.RemoveField(
            model_name="paper",
            name="citations",
        ),
        migrations.RemoveField(
            model_name="paper",
            name="influential_citations",
        ),
    ]
