# Generated by Django 5.0.1 on 2024-01-25 11:07

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("vetu", "0011_authorimpact_last_updated"),
    ]

    operations = [
        migrations.AddField(
            model_name="paper",
            name="topic_codes",
            field=models.CharField(
                blank=True, default="Not Coded Yet", max_length=100, null=True
            ),
        ),
    ]