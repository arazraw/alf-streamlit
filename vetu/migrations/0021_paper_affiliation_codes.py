# Generated by Django 5.0.4 on 2024-04-23 06:34

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("vetu", "0020_alter_paper_journal_title_and_more"),
    ]

    operations = [
        migrations.AddField(
            model_name="paper",
            name="affiliation_codes",
            field=models.TextField(default=None),
            preserve_default=False,
        ),
    ]