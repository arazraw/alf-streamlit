# Generated by Django 4.2.9 on 2024-02-05 10:38

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("vetu", "0016_alter_author_identifier"),
    ]

    operations = [
        migrations.AlterField(
            model_name="author",
            name="name",
            field=models.CharField(max_length=200),
        ),
    ]
