# Generated by Django 5.0.1 on 2024-01-17 07:41

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("vetu", "0001_initial"),
    ]

    operations = [
        migrations.AddField(
            model_name="paper",
            name="affiliations",
            field=models.TextField(default=False),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name="paper",
            name="affiliations_checked",
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name="paper",
            name="akademisk",
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name="paper",
            name="got_un",
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name="paper",
            name="kar_in",
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name="paper",
            name="lin_un",
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name="paper",
            name="lun_un",
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name="paper",
            name="ore_un",
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name="paper",
            name="region",
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name="paper",
            name="ume_un",
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name="paper",
            name="upp_un",
            field=models.BooleanField(default=False),
        ),
    ]
