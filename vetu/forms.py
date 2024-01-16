from django import forms
from vetu.models import Paper

class PaperFilterForm(forms.Form):
    start_year = forms.IntegerField(required=False, min_value=1990, max_value=2024, label='Start Year')
    end_year = forms.IntegerField(required=False, min_value=1990, max_value=2024, label='End Year')