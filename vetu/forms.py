from django import forms
from vetu.models import Paper

class PaperFilterForm(forms.Form):
    start_year = forms.IntegerField(required=False)
    end_year = forms.IntegerField(required=False)
    akademisk = forms.BooleanField(required=False)
    region = forms.BooleanField(required=False)