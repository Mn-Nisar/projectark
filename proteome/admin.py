from django.contrib import admin
from . models import DataAnalysis , Example , Contaminant , FastaBackup , Ploting , Forzip, PlotExample

@admin.register(DataAnalysis)
class DataAnalysisAdmin(admin.ModelAdmin):
    list_display = ['user']


@admin.register(Example)
class ExampleAdmin(admin.ModelAdmin):
    pass

@admin.register(Contaminant)
class ContaminantAdmin(admin.ModelAdmin):
    pass

@admin.register(FastaBackup)
class FastaBackupAdmin(admin.ModelAdmin):
    pass

@admin.register(Ploting)
class PlotingBackupAdmin(admin.ModelAdmin):
    pass

@admin.register(Forzip)
class ForzipAdmin(admin.ModelAdmin):
    pass


@admin.register(PlotExample)
class ForzipAdmin(admin.ModelAdmin):
    pass
