# Generated by Django 4.2.9 on 2024-02-10 09:27

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("vetu", "0018_alter_paper_topic"),
    ]

    operations = [
        migrations.AlterField(
            model_name="paper",
            name="topic_codes",
            field=models.CharField(
                blank=True, default="Not Coded Yet", max_length=600, null=True
            ),
        ),
    ]