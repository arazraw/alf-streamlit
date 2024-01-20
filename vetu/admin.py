from django.contrib import admin
from .models import Author, Paper, Impact, AuthorImpact

admin.site.register(Author)
admin.site.register(Paper)
class ImpactAdmin(admin.ModelAdmin):
    list_display = ('paper', 'citations', 'impactful_citations', 'meaningful', 'last_updated')

admin.site.register(Impact, ImpactAdmin)
admin.site.register(AuthorImpact)