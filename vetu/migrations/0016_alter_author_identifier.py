# Generated by Django 4.2.9 on 2024-02-05 10:35

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("vetu", "0015_alter_author_affiliation"),
    ]

    operations = [
        migrations.AlterField(
            model_name="author",
            name="identifier",
            field=models.CharField(blank=True, max_length=500, null=True),
        ),
    ]
