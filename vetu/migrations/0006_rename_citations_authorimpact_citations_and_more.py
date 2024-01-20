# Generated by Django 5.0.1 on 2024-01-19 09:31

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("vetu", "0005_authorimpact"),
    ]

    operations = [
        migrations.RenameField(
            model_name="authorimpact",
            old_name="Citations",
            new_name="citations",
        ),
        migrations.RenameField(
            model_name="authorimpact",
            old_name="Impactful_Citations",
            new_name="impactful_Citations",
        ),
        migrations.RenameField(
            model_name="authorimpact",
            old_name="Name",
            new_name="name",
        ),
        migrations.RenameField(
            model_name="authorimpact",
            old_name="ORCID",
            new_name="orcid",
        ),
        migrations.AddField(
            model_name="authorimpact",
            name="paper_count",
            field=models.IntegerField(default=0),
        ),
        migrations.AddField(
            model_name="authorimpact",
            name="s2_id",
            field=models.CharField(default=0, max_length=20),
            preserve_default=False,
        ),
    ]