# Generated by Django 5.0.1 on 2024-01-15 20:49

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("search", "0008_impact_last_updated"),
    ]

    operations = [
        migrations.AlterField(
            model_name="impact",
            name="citation_data",
            field=models.JSONField(default=dict, null=True),
        ),
    ]
