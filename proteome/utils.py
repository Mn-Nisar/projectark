from django.core.files.images import ImageFile
from django.core.files.base import ContentFile

import pandas as pd
from difflib import SequenceMatcher
import re
from itertools import chain
import base64 
import io
import stringdb
# from my_gprofiler import GProfiler
from . import my_gprofiler
from .models import DataAnalysis , Ploting , PlotExample
from .species import ALL_SPECIES

import time
import multiprocessing

def get_protien_interaction(genes):
    string_ids = stringdb.get_string_ids(genes)
    enrichment_df = stringdb.get_interaction_partners(string_ids.queryItem , required_score=900,)
    return enrichment_df

def abundances(columns):
    abundance_list =  []
    for l in columns:
        l = l.strip()
        if ('Abundances' in l) or ('Abundance' in l):
            abundance_list.append(l)

    final_list = list()
    for abd in abundance_list:
        if ("Ratio" in abd) or ("Grouped" in abd) or ("Normalized" in abd) or ("normalized" in abd):
            continue
        else:
            final_list.append(abd)
    
    if  len(final_list) > 0:  
        final_list.sort()
        return final_list
    else:
        return columns

    #give option to see all column names also if the name column name may vary

def clean_coulumn_heading(sample_data_columns):

    samples = sample_data_columns.split("SaMpSepeR")

    sample_list = []
    for sample in samples:
        abd_list = []
        abundance = sample.split('RepsepRatTor')
        for abd in abundance:
            abd = abd.strip()   
            if abd != '':
                abd_list.append(abd)
        sample_list.append(abd_list)
    final_list = [x for x in sample_list if x != []]
    return final_list

def sort_name(samples):
    new_samp = []
    for s in samples:
        if 'Abundances' in s:
            s = s.replace('Abundances','')
        if  'Abundance' in s:
            s = s.replace('Abundance','')
        if  'NORM_' in s:
            s = s.replace('NORM_','')
        if 'AVG_NORM' in s:
            s = s.replace('AVG_NORM','')
        s = s.strip()
        new_samp.append(s)


    name = ''
    if len(new_samp) > 1:
        string1 = new_samp[0]
        string2 = new_samp[1]
        match = SequenceMatcher(None, string1,
                        string2).find_longest_match(0, len(string1), 0,
                                                    len(string2))
        name =  string1[match.a:match.a + match.size]
    else:
        name = samples

    return name


def removeSpaceAndComma(columns):
    cleaned_col = []
    for column in columns:
            if ',' in column:
                column =column.strip()
                column = column.replace(',',' ')
                cleaned_col.append(column)
            else:
                column =column.strip()
                cleaned_col.append(column)

    return cleaned_col

def forPCA(sample_columns,control_columns,sna,cna):
    before_norm = []
    after_norm = []
    for sample_list in sample_columns:
        for sample in sample_list:
            before_norm.append(sample)

    for control_list in control_columns:
        for control in control_list:
            before_norm.append(control)

    for norm_sample_list in sna:
        for sample in norm_sample_list:
            after_norm.append(sample)

    for norm_control_list in cna:
        for control in norm_control_list:
            after_norm.append(control)

    return before_norm,after_norm


def expandNCleanColumns(sample_columns,control_columns):
    colum_list_sample = []
    for sample_list in sample_columns:
        for sample in sample_list:
            colum_list_sample.append(sample)

    column_list_control = []
    for control_list in control_columns:
        for control in control_list:
            column_list_control.append(control)

    colum_list_sample = removeSpaceAndComma(colum_list_sample)
    column_list_control = removeSpaceAndComma(column_list_control)

    return colum_list_sample, column_list_control

def expandcols(sample_columns,control_columns):

    sample = list(chain.from_iterable(sample_columns))
    control = list(chain.from_iterable(control_columns))
    sna = ["NORM_"+x for x in sample]
    cna = ["NORM_"+x for x in control]
    return sample, control, sna , cna


def expandforQuant(sample_columns,control_columns):
    colum_list_sample = []
    sna = []
    for sample_list in sample_columns:
        normlist = []
        for sample in sample_list:
            colum_list_sample.append(sample)
            normlist.append("NORM_"+sample)
        sna.append(normlist)

    cna = []
    column_list_control = []
    
    for control_list in control_columns:
        contnormlist = []
        for control in control_list:
            column_list_control.append(control)
            contnormlist.append("NORM_"+control)
        cna.append(contnormlist)
    return colum_list_sample,column_list_control,sna, cna


def intensities(columns):
    intensitiy_list =  []
    for l in columns:
        l = l.strip()
        if ('intensity' in l) or ('Intensity' in l) or ('intensities' in l) or ('Intensities' in l) or ('Abundances' in l) or ('Abundance' in l):
            intensitiy_list.append(l)
    intensitiy_list.sort()
    return intensitiy_list

def removeavgsmp(avg_sample):
    if 'AVG_NORM' in avg_sample:
        colname = avg_sample.replace('AVG_NORM','')
        return colname
    else:
        return avg_sample

def lablesforbox(columns):
    labels = {}
    for column in columns:
         labels[column] = truncc(column)
    return labels

def truncc(column):
    column = column.replace('Abundance','')
    column = column.replace('Abundances','')
    column = column.replace('NORM_','')
    column = column.split()
    return column[-1]

def columnsforbox(columns):

    clean_col = list()
    for column in columns:
        if 'Abundances' in column:
            column = column.replace('Abundances','')
        if 'Abundance' in column:
            column = column.replace('Abundance','')
        if 'NORM_' in column:
            column = column.replace('NORM_','')

        clean_col.append(column)

    return clean_col



def getAccesion(columns):
    keys = []
    for cols in columns:
        accession = re.search('.*[aA][cC][cC][eE][sS][sS][iI][oO][nN].*',cols)
        if accession:
            keys.append(cols)
        accession = re.search('.*[pP][rR][oO][tT][eE][in][nN].*',cols)
        if accession:
            keys.append(cols)
    return keys



def seperatesamples(sample_columns,control_columns):

    sortedsample = []
    for sample in sample_columns:
        sample.sort()
        sortedsample.append(sample)

    i = len(sortedsample[0])

    final = list()
    sna = list()
    batchlist = list()

    for j in range(i):
        newsamp = []
        newsampn = []
        for key, sample in  enumerate(sortedsample):
            newsamp.append(sample[j])
            newsampn.append("NORM_"+sample[j])
            batchlist.append(key+1)
        final.append(newsamp)
        sna.append(newsampn)


    sortedcontrol = []
    for sample in control_columns:
        sample.sort()
        sortedcontrol.append(sample)

    i = len(sortedcontrol[0])

    final_control = list()
    cna = list()
    batchlist = list()

    for j in range(i):
        newcont = []
        newcontn = []
        for key, sample in  enumerate(sortedcontrol):
            newcont.append(sample[j])
            newcontn.append("NORM_"+sample[j])
            batchlist.append(key+1)
        final_control.append(newcont)
        cna.append(newcontn)

    return final, final_control, sna, cna


def clean_custom_names(sample):
    cleaned_s = list()
    nameDict = dict()
    samples = sample.split('ZohaNSP')
    for s in samples:
        s = s.strip()
        cleaned_s.append(s)
    cleaned_s = [x for x in cleaned_s if x]

    for n, s in enumerate(cleaned_s):
        nameDict[n] = s
    return nameDict


def gene_ontology_calc(accession,species, pvcutOff ):

    selected_species = ALL_SPECIES[species]
    try:
        gp = my_gprofiler.GProfiler(return_dataframe=True)
        gene_ont = gp.profile(organism=selected_species,query=accession,no_evidences=False, user_threshold = 1)
    except:
        gp = my_gprofiler.GProfiler(return_dataframe=True)
        gene_ont = gp.profile(organism=selected_species,query=accession,no_evidences=False, user_thresh = pvcutOff)

    gene_ont.sort_values(by=['source'], inplace=True)
    return gene_ont


def cobmine_samp_control(samp_col, control_col):
    new_columns = list()
    i = 0
    for sample in samp_col:
        each_btch = sample + control_col[i]
        new_columns.append(each_btch)
        i += 1

    fordf = list(chain.from_iterable(new_columns))
    return new_columns, fordf


def clean_plot_df(columns):
    newcols = [x.replace(',',' ').strip() for x in columns]
    return newcols


def get_gene_symbol(columns):

    keys = []
    for cols in columns:
        gene = re.search('.*[gG][eE][nN][eE].*',cols)
        if gene:
            keys.append(cols)
    return keys


def convert_acc_to_gene(accessions):
    
    accessions = [acc.split(';')[0] if ";" in acc else acc for acc in accessions]
    
    gs_convert_success = True
    con_df = pd.DataFrame()
    try:
        gp = my_gprofiler.GProfiler(return_dataframe=True)
        
        gs = gp.convert(organism='hsapiens',
                query=accessions,
                target_namespace='ENTREZGENE_ACC')

        con_df = gs[['incoming','name']]
        con_df.rename(columns = {'incoming':'Accesion_gf' , 'name': '_GENE_SYMBOL_'}, inplace = True)
        
        none_df = con_df.loc[ con_df['_GENE_SYMBOL_'] == 'None']

        for index, row in none_df.iterrows():
            if row['_GENE_SYMBOL_'] == 'None':
                con_df.loc[index, '_GENE_SYMBOL_'] =  con_df.loc[index, 'Accesion_gf']

    except:
        gs_convert_success = False

    return con_df , gs_convert_success


def calc_avg_exp(x):

    x = x.dropna()
    x = list(x)
    upcount = x.count("Upregulated")
    downcount = x.count("Downregulated")

    if upcount >= downcount:
        return "Upregulated"
    else:
        return "Downregulated"

def stringtolist(x):
    x = x.strip("][")
    x = x.split(',')
    listt = [i.replace("'",'').strip() for i in x]
    return listt


def decode_svg_for_zip(file_code , name):
    decoded = base64.b64decode(file_code)
    decoded_image  = ImageFile(io.BytesIO(decoded), name = name)
    return decoded_image




def savefile_csv(files , user ,lablled , lablefree, diffex):

    if files.name.endswith('.xlsx'):
        df = pd.read_excel(files,engine='openpyxl')
        df.columns = clean_plot_df(list(df.columns))
        columns  = df.columns
        df = df.to_csv(index = False)
        updated_file = ContentFile(df)
        updated_file.name = "inputfile.csv"
        
        if diffex:
            data_als = DataAnalysis.objects.create(resultData = updated_file, user = user, labledData = lablled,lableFreeData = lablefree )
        else:
            data_als = DataAnalysis.objects.create(file = updated_file, user = user, labledData = lablled,lableFreeData = lablefree )
        
        job_id = data_als.id
        data_als.save()
        


    elif files.name.endswith('.xls'):
        df = pd.read_excel(files,engine='xlrd')
        df.columns = clean_plot_df(list(df.columns))

        columns  = df.columns

        df = df.to_csv(index = False)
        updated_file = ContentFile(df)
        updated_file.name = "inputfile.csv"

        if diffex:
            data_als = DataAnalysis.objects.create(resultData = updated_file, user = user, labledData = lablled,lableFreeData = lablefree )

        else:
            data_als = DataAnalysis.objects.create(file = updated_file, user = user, labledData = lablled,lableFreeData = lablefree )
   

        job_id = data_als.id
        data_als.save()
        


    elif files.name.endswith('.csv'):

        df = pd.read_csv(io.StringIO(files.read().decode('utf-8')) )
        df.columns = clean_plot_df(list(df.columns))
        columns  = df.columns
        df = df.to_csv(index = False)
        updated_file = ContentFile(df)
        updated_file.name = "inputfile.csv"
    
        if diffex:
            data_als = DataAnalysis.objects.create(resultData = updated_file, user = user, labledData = lablled,lableFreeData = lablefree )

        else:
            data_als = DataAnalysis.objects.create(file = updated_file, user = user, labledData = lablled,lableFreeData = lablefree )
   
        job_id = data_als.id
        data_als.save()
        
    elif files.name.endswith('.txt') or files.name.endswith('.tsv'):

        df = pd.read_csv(io.StringIO(files.read().decode('utf-8')) , delimiter='\t')
        df.columns = clean_plot_df(list(df.columns))
        columns  = df.columns
        df = df.to_csv(index = False)
        updated_file = ContentFile(df)
        updated_file.name = "inputfile.csv"

        if diffex:
            data_als = DataAnalysis.objects.create(resultData = updated_file, user = user, labledData = lablled,lableFreeData = lablefree )

        else:
            data_als = DataAnalysis.objects.create(file = updated_file, user = user, labledData = lablled,lableFreeData = lablefree )
           
        job_id = data_als.id
        data_als.save()

    return job_id , list(columns)


def getAccesionCol(columns):

    keys = []
    for cols in columns:
        accession = re.search('.*[aA][cC][cC][eE][sS][sS][iI][oO][nN].*',cols)
        if accession:
            keys.append(cols)
        accession = re.search('.*[pP][rR][oO][tT][eE][in][nN].*',cols)
        if accession:
            keys.append(cols)

        accession = re.search('.*[iI][dD].*',cols)
        if accession:
            keys.append(cols)

        accession = re.search('.*[gG][eE][nN][eE].*',cols)
        if accession:
            keys.append(cols)
    
    final = keys + [x for x in columns if x not in keys]

    return final


def getGeneCol(columns):

    keys = []
    for cols in columns:
        accession = re.search('.*[gG][eE][nN][eE].*',cols)
        if accession:
            keys.append(cols)

        accession = re.search('.*[iI][dD].*',cols)
        if accession:
            keys.append(cols)

        accession = re.search('.*[aA][cC][cC][eE][sS][sS][iI][oO][nN].*',cols)
        if accession:
            keys.append(cols)
        accession = re.search('.*[pP][rR][oO][tT][eE][in][nN].*',cols)
        if accession:
            keys.append(cols)

    final = keys + [x for x in columns if x not in keys]

    return final




def save_plot_file(files, plot_type , user):

    if files.name.endswith('.xlsx'):
        df = pd.read_excel(files,engine='openpyxl')
        df.columns = clean_plot_df(list(df.columns))
        columns_temp  =  df.columns

        column_dtype = df.dtypes.to_dict()
        df = df.to_csv(index = False)
        updated_file = ContentFile(df)
        updated_file.name = "inputfile.csv"
        plot_qs = Ploting.objects.create(file = updated_file,plot_type = plot_type)
        plot_id = plot_qs.id
        plot_qs.save()

        return plot_id , columns_temp , column_dtype

    elif files.name.endswith('.xls'):
        df = pd.read_excel(files,engine='xlrd')
        df.columns = clean_plot_df(list(df.columns))
        columns_temp  =  df.columns
        column_dtype = df.dtypes.to_dict()

        df = df.to_csv(index = False)
        updated_file = ContentFile(df)
        updated_file.name = "inputfile.csv"
        plot_qs = Ploting.objects.create(file = updated_file , plot_type = plot_type)
        plot_id = plot_qs.id
        plot_qs.save()
        

        return plot_id , columns_temp , column_dtype


    elif files.name.endswith('.csv'):
        
        df = pd.read_csv(io.StringIO(files.read().decode('utf-8')) )
    
        df.columns = clean_plot_df(list(df.columns))

        column_dtype = df.dtypes.to_dict()

        columns_temp  =  df.columns
        df = df.to_csv(index = False)
        updated_file = ContentFile(df)
        updated_file.name = "inputfile.csv"
        plot_qs = Ploting.objects.create(file = updated_file , plot_type = plot_type)
        plot_id = plot_qs.id
        plot_qs.save()
        
        return plot_id  , columns_temp , column_dtype

    elif files.name.endswith('.txt') or files.name.endswith('.tsv'):

        df = pd.read_csv(io.StringIO(files.read().decode('utf-8')) , delimiter='\t')
    
        df.columns = clean_plot_df(list(df.columns))
    
        column_dtype = df.dtypes.to_dict()

        columns_temp  =  df.columns
        df = df.to_csv(index = False)
        updated_file = ContentFile(df)
        updated_file.name = "inputfile.csv"
        plot_qs = Ploting.objects.create(file = updated_file , plot_type = plot_type)
        plot_id = plot_qs.id
        plot_qs.save()

        return plot_id  , columns_temp , column_dtype

    else:
        return None

def load_explot_file(plot_type , user):
    
    q = PlotExample.objects.get(plot_type = plot_type)
    
    files = q.file.path
    
    if files.endswith('.xlsx'):
        df = pd.read_excel(files,engine='openpyxl')
        df.columns = clean_plot_df(list(df.columns))
        columns_temp  =  df.columns
        
        column_dtype = df.dtypes.to_dict()

        df = df.to_csv(index = False)
        updated_file = ContentFile(df)
        updated_file.name = "inputfile.csv"
        plot_qs = Ploting.objects.create(file = updated_file , plot_type = plot_type)
        plot_id = plot_qs.id
        plot_qs.save()
        
        return plot_id , columns_temp , column_dtype

    elif files.endswith('.xls'):
        df = pd.read_excel(files,engine='xlrd')

        df.columns = clean_plot_df(list(df.columns))

        columns_temp  =  df.columns
        
        column_dtype = df.dtypes.to_dict()

        df = df.to_csv(index = False)
        updated_file = ContentFile(df)
        updated_file.name = "inputfile.csv"
        plot_qs = Ploting.objects.create(file = updated_file, user = user , plot_type = plot_type)
        plot_id = plot_qs.id
        plot_qs.save()
        
        return plot_id , columns_temp , column_dtype


    elif files.endswith('.csv'):
        
        df = pd.read_csv(files)
    
        df.columns = clean_plot_df(list(df.columns))
        columns_temp  =  df.columns

        column_dtype = df.dtypes.to_dict()

        df = df.to_csv(index = False)
        updated_file = ContentFile(df)
        updated_file.name = "inputfile.csv"
        plot_qs = Ploting.objects.create(file = updated_file, user = user , plot_type = plot_type)
        plot_id = plot_qs.id
        plot_qs.save()
        
        return plot_id  , columns_temp , column_dtype

    elif files.endswith('.txt') or files.endswith('.tsv'):

        df = pd.read_csv( files , delimiter='\t')
    
        df.columns = clean_plot_df(list(df.columns))
        columns_temp  =  df.columns
        
        column_dtype = df.dtypes.to_dict()

        df = df.to_csv(index = False)
        updated_file = ContentFile(df)
        updated_file.name = "inputfile.csv"
        plot_qs = Ploting.objects.create(file = updated_file, user = user , plot_type = plot_type)
        plot_id = plot_qs.id
        plot_qs.save()
        
        return plot_id  , columns_temp , column_dtype

    else:
        return None


def get_sna_cna(samp_col , con_col):
    sna = [["NORM_"+cols for cols in i] for i in samp_col]
    cna = [["NORM_"+cols for cols in i] for i in con_col]
    return sna , cna 

def  getbatches(s_col , c_col):
    batch_list = []
    sortedsample = []
    sortedcontrol = []
    final_cols = []

    for sample in s_col:
        sample.sort()
        sortedsample.append(sample)

    for control in c_col:
        control.sort()
        sortedcontrol.append(control)

    i = 1
    for batch in sortedsample:
        for sample in batch:
            batch_list.append(i)
            final_cols.append(sample)
        i +=1

    i = 1
    for batch in sortedcontrol:
        for control in batch:
            batch_list.append(i)
            final_cols.append(control)
        i +=1

    return batch_list, final_cols



def get_matrix_limma(sna,controlcols,snames):    
    matrix = []
    fordf = []
    pval_array = []
    fc_array = []
    rename_list = []
    
    matrix.append(["control" for c in range(0,len(controlcols))])
    fordf.append(controlcols)

    for key,s in snames.items():
        
        pval_array.append('P VALUE of'+s)
        fc_array.append('FOLDCHANGE_'+s)

        s_for_matrix = ['sample'+key for i in range(0,len( sna[int(key)] ))]
        matrix.append(s_for_matrix)

    matrix = list(chain.from_iterable(matrix))

    fordf.append(chain.from_iterable(sna))
    fordf = list(chain.from_iterable(fordf))

    for k,v in enumerate(pval_array):
        rename_list.append(v)
        rename_list.append(fc_array[k])

    return matrix,fordf,pval_array, fc_array , rename_list
