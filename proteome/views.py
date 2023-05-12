from django.shortcuts import render 
from django.http import HttpResponse
from django.http import Http404
from django.core.files.base import ContentFile

from django.conf import settings
from django.contrib import messages
from django.http import JsonResponse

import io
import os
import zipfile

import pandas as pd
import numpy as np
import matplotlib


from .utils import (abundances, savefile_csv , clean_coulumn_heading,intensities,lablesforbox, clean_custom_names,
    sort_name , columnsforbox, seperatesamples,cobmine_samp_control, calc_avg_exp, gene_ontology_calc, getAccesion,get_protien_interaction, clean_plot_df, get_gene_symbol,
     stringtolist, decode_svg_for_zip, getAccesionCol, getGeneCol, convert_acc_to_gene , save_plot_file, load_explot_file,  getbatches)

from .keggpath import draw_pathway
from .models import DataAnalysis , Example, Contaminant , Ploting , Forzip, PlotExample

def home(request):
    matplotlib.rc_file_defaults()
    return render(request,'proteome/index.html')

def input(request):
    return render(request,'proteome/home.html')

def inputf(request): 
    lablled = True
    lablefree = False

    if (request.method == 'POST'):
        example_analysis = request.POST.get("wokring-with-ex")
        
        print(example_analysis)

        if example_analysis == "yes":
            
            q = Example.objects.filter(usethis = True).first()
            files = q.file
        else:
            files = request.FILES['file']

        
        
        if files.name.endswith('.xlsx') or files.name.endswith('.csv') or files.name.endswith('.txt') or files.name.endswith('.xls'):

            user = None
            job_id , columns = savefile_csv(files , user ,lablled , lablefree, False)
            
            request.session['job_id'] = job_id

            accession_col = getAccesionCol(columns)
        
            gene_col = getGeneCol(columns)

            if (request.POST.get('rep_method')) == "techrep":

                number_of_samples = int(request.POST.get('no_of_sample'))
                number_of_control = int(request.POST.get('no_of_control'))

                abd_columns = abundances(columns)
                
                context = {'abd_columns':abd_columns,'number_of_samples':number_of_samples,
                'number_of_control':number_of_control, 'accession_col':accession_col, 'gene_col':gene_col}
                
                return render(request,'proteome/pre_analyze.html',context)

            else:
                number_of_batches = int(request.POST.get('no_of_batches'))
                samples_in_bio = request.POST.get('samples_in_bio')
                request.session['samples_in_bio'] = samples_in_bio

                abd_columns = abundances(columns)

                context = {'abd_columns':abd_columns,'number_of_batches':number_of_batches,
                'number_of_samples':samples_in_bio, 'accession_col':accession_col, 'gene_col':gene_col}

                return render(request,'proteome/pre_anlz_bio.html',context)

        else:
            messages.error(request, 'Please upload only Excel , CSV or TXT file')
            return render(request, 'proteome/home.html')

    return render(request, 'proteome/home.html')