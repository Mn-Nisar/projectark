from django.urls import path
from django.views.generic import TemplateView
from . import views

app_name = 'proteome'

urlpatterns = [
    path('',views.home, name = 'home'),
    path('input/', views.input, name = 'input'),
    path('inputfile/', views.inputf, name='inputfile'),

    
]


