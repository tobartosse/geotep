#!/usr/bin/env python3

def ComparativeMatrix():
    files =  os.listdir('./ExpPro')
    for f in files:
        if '_fc.tsv' in f:
            name = f.split('_fc')[0]
            out1 = open(f'./ExpPro/{name}_Comparative_Matrix_FC.txt', 'w')
            file = open(f'./ExpPro/{f}')
            rfile = file.readline()
            out1.write(f'{rfile}')
            rfile = file.readline()
            while rfile!='':
                line = rfile.strip().split()
                check = line[1:]
                pss = 'yes'
                for i in check:
                    i=float(i)
                    if i==0:
                        pss='no'
                if pss == 'yes':
                    out1.write(f'{rfile}')

                rfile = file.readline()

        if '_tt.tsv' in f:
            name = f.split('_tt')[0]
            out2 = open(f'./ExpPro/{name}_Comparative_Matrix_TT.txt', 'w')
            file = open(f'./ExpPro/{f}')
            rfile = file.readline()
            out2.write(f'{rfile}')
            rfile = file.readline()
            while rfile!='':
                line = rfile.strip().split()
                check = line[1:]
                pss = 'yes'
                for i in check:
                    i=float(i)
                    if i==0:
                        pss='no'
                if pss == 'yes':
                    out2.write(f'{rfile}')
                rfile = file.readline()

def Pandas(f, t, e, dicinfo,lg):
    df = pd.read_csv(f'./ALL/{f}', sep = '\t', index_col='ID_REF')
    for k in dicinfo[t][e]:
        C = dicinfo[t][e][k]['C']
        S = dicinfo[t][e][k]['S']
        #Calculate Log2FC
        df['aveC'] = df[C].mean(axis=1)
        df['aveS'] = df[S].mean(axis=1)
        df = df[(df != 0).all(axis=1)]####
        df['Log2FC'] = np.log2(df['aveS']/df['aveC'])
        #Application to T-TEST
        def apply_t_test(row):
            samples = row[S]
            controls = row[C]
            t_statistic, p_value = stats.ttest_ind(samples, controls, alternative='greater')
            return p_value
        df['p_value'] = df.apply(apply_t_test, axis=1)
    
        target_gene = lg#OJO ESTE SE DEFINE DESDE EL INICIO EN UNA LISTA DE GENES
        TG_FC = df.at[target_gene, 'Log2FC']
        if TG_FC>0:
            df_k = df[(df['Log2FC']>0)]
            df_std = df_k['Log2FC'].std()
        if TG_FC<0:
            df_k = df[(df['Log2FC']<0)]
            df_std = df_k['Log2FC'].std()
        inf = TG_FC-(df_std/10)
        sup = TG_FC+(df_std/10)
        #print (inf, TG_FC, sup)

        df_TG = df[(df['Log2FC'] > inf) & (df['Log2FC'] <= sup)]
        #to get only the calcs
        calcs = ['Log2FC', 'p_value']
        df_TG = df_TG[calcs]
        
        df_TG.to_csv(f'./temp/{t}_{e}_{k}.TG.FCTT.csv', index=True)
        
def getdata(lg):
    file = open('./ALL/Info_samples.txt')
    rfile =  file.readline()
    dicinfo = {}
    while rfile!='':
        l = rfile.strip().split()
        t = l[0]
        ex = l[1]
        sa = l[2]
        cl = l[3:]
        if t not in dicinfo:
            dicinfo[t] = {}
            if ex not in dicinfo[t]:
                dicinfo[t][ex] = {}
                if sa not in dicinfo[t][ex]:
                    dicinfo[t][ex][sa] = {}
                    for s in cl:
                        sname = s.split('_')[0]
                        stype = s.split('_')[1]
                        if stype not in dicinfo[t][ex][sa]:
                            dicinfo[t][ex][sa][stype] = [sname]
                        else:
                            dicinfo[t][ex][sa][stype] += [sname]
            else:
                if sa not in dicinfo[t][ex]:
                    dicinfo[t][ex][sa] = {}
                    for s in cl:
                        sname = s.split('_')[0]
                        stype = s.split('_')[1]
                        if stype not in dicinfo[t][ex][sa]:
                            dicinfo[t][ex][sa][stype] = [sname]
                        else:
                            dicinfo[t][ex][sa][stype] += [sname]
        else:
            if ex not in dicinfo[t]:
                dicinfo[t][ex] = {}
                if sa not in dicinfo[t][ex]:
                    dicinfo[t][ex][sa] = {}
                    for s in cl:
                        sname = s.split('_')[0]
                        stype = s.split('_')[1]
                        if stype not in dicinfo[t][ex][sa]:
                            dicinfo[t][ex][sa][stype] = [sname]
                        else:
                            dicinfo[t][ex][sa][stype] += [sname]
            else:
                if sa not in dicinfo[t][ex]:
                    dicinfo[t][ex][sa] = {}
                    for s in cl:
                        sname = s.split('_')[0]
                        stype = s.split('_')[1]
                        if stype not in dicinfo[t][ex][sa]:
                            dicinfo[t][ex][sa][stype] = [sname]
                        else:
                            dicinfo[t][ex][sa][stype] += [sname]
        rfile = file.readline()

    #for t in dicinfo:
    #    for e in dicinfo[t]:
    #        for s in dicinfo[t][e]:
    #            for v in dicinfo[t][e][s]:               
    #                print (t, e, s, v, dicinfo[t][e][s][v])

    #load files
    files = os.listdir('./ALL')
    for f in files:
        for t in dicinfo:
            for e in dicinfo[t]:
                if e in f:
                    Pandas(f, t, e, dicinfo, lg)    

def makematrix(lg):
    #load name
    file = open('./GEO/genename.txt')
    rfile = file.readline()
    dicname = {}
    while rfile!='':
        n1 = rfile.strip().split()[0]
        n2 = rfile.strip().split()[1]
        dicname[n1]=n2
        rfile = file.readline()

    files = os.listdir('./temp')
    dicFC = {}
    dicTT = {}
    dicpst, n = {}, 0
    for f in files:
        name = f.strip().split('.')[0]
        dicpst[name]=n
        n+=1
        file = open(f'./temp/{f}')
        rfile = file.readline()
        rfile = file.readline()
        while rfile!='':
            l=rfile.strip().split(',')
            gene = l[0]
            ngene = gene
            if gene in dicname:
                ngene=dicname[gene]
            fc = float(l[1])
            tt = float(l[2])
            if ngene not in dicFC:
                dicFC[ngene]=[0]*len(files)
                dicTT[ngene]=[0]*len(files)
                dicFC[ngene][dicpst[name]]=fc
                dicTT[ngene][dicpst[name]]=tt
            else:
                dicFC[ngene][dicpst[name]]=fc
                dicTT[ngene][dicpst[name]]=tt

            rfile = file.readline()
    invdic = {}
    for f in dicpst:
        invdic[dicpst[f]]=f
    labels = 'Gene'
    for i in range(len(dicpst)):
        labels+=f'\t{invdic[i]}'
    outfc = open(f'./ExpPro/{lg}_fc.tsv', 'w')
    outtt = open(f'./ExpPro/{lg}_tt.tsv', 'w')
    outfc.write(f'{labels}\n')
    outtt.write(f'{labels}\n')
    
    for g in dicFC:
        values = ''
        for v in dicFC[g]:
            values+=f'\t{v}'
        outfc.write(f'{g}{values}\n')
    outfc.close()

    for g in dicTT:
        values = ''
        for v in dicTT[g]:
            values+=f'\t{v}'
        outtt.write(f'{g}\t{values}\n')
    outtt.close()
def analysis(addfile):
    df = pd.read_csv(f'{addfile}', sep = '\t', skiprows=3)
    df = df.rename(columns={'!UserClassification': 'GENE'})
    df.set_index('GENE', inplace=True)
    
    #Get Samples and Controls
    C = [col for col in df.columns if 'C' in col]
    S = [col for col in df.columns if 'S' in col]
    N = [col for col in df.columns if 'N' in col]
  
    #Calculate Log2FC
    df['aveC'] = df[C].mean(axis=1)
    df['aveS'] = df[S].mean(axis=1)
    df['Log2FC'] = np.log2(df['aveS']/df['aveC'])

    #Application to T-TEST
    def apply_t_test(row):
        samples = row[S]
        controls = row[C]
        t_statistic, p_value = stats.ttest_ind(samples, controls, alternative='greater')
        return p_value
    
    df['p_value'] = df.apply(apply_t_test, axis=1)

    #print (df)
    
    #GET MATRIX FOR ONLY GENES WITH THE SAME VALUES OF A TARGET GENES
    outname = addfile[:-7]
    #target_gene = input(f"Check the soft file for {outname}, and write here the name of the target gene, ex:    : ")
    target_gene = '262234_at'
    TG_FC = df.at[target_gene, 'Log2FC']

    if TG_FC>0:
        df_k = df[(df['Log2FC']>0)]
        df_std = df_k['Log2FC'].std()
    if TG_FC<0:
        df_k = df[(df['Log2FC']<0)]
        df_std = df_k['Log2FC'].std()
    inf = TG_FC-(df_std/10)
    sup = TG_FC+(df_std/10)
    #print (inf, TG_FC, sup)

    df_TG = df[(df['Log2FC'] > inf) & (df['Log2FC'] <= sup)]
    #to get only the calcs
    calcs = ['Log2FC', 'p_value']
    df_TG = df_TG[calcs]

    df_TG.to_excel(f'{outname}.TG.FCTT.xlsx')

    ##OJO AQUI SE GUARDAN LOS RESULTADOS
    #Qui apago las siguientes lineas porque son cuando se corre por defecto
    #Como estamos buscando un gen especifico, esto archivos no los necesito!
    #pero OJO, si estas corriendo otrso experimientos Aplica completamente!
    """
    HEG = df[(df['Log2FC'] > 1.5) & (df['p_value'] < 0.05)]
    LEG = df[(df['Log2FC'] < -1.5) & (df['p_value'] < 0.05)]

    
    HEG.to_excel(f'{outname}.HEG.xlsx')
    HEG.to_excel(f'{outname}.LEG.xlsx')
    """

def loadmatrix():
    out = open('./ALL/Info_samples.txt', 'w')
    folders = os.listdir('./GEO')
    for folder in folders:
        if '.' not in folder:
            files = os.listdir(f'./GEO/{folder}')
            for file in files:
                #make review samples
                if 'mf.tsv' in file:
                    exp = file.split('.')[0]
                    info = input (f"""
                    IMPORTANT: 
                    Give a simple name for the experiment {exp} found in the folder {folder}, 
                    it allows to identify the kind of comparison and profile (for the next step)
                    that is calculated by using these comparable scores. Write the name,
                    Ex. Tolerant, Resistant or Mutant, etc : """)
                    print (f"""
                    Calculating Log2FC and T-Test scores for {exp} experiment or {info} comparison...
                    """)
                    f = open(f'./GEO/{folder}/{file}')
                    rf = f.readlines()
                    id = rf[2].strip().split('\t')[1:]
                    tid = rf[3].strip().split('\t')[1:]
                    idfinal = ''
                    for i in range(len(id)):
                        idfinal+=f'\t{id[i][1:-1]}_{tid[i]}'
                    out.write(f'{folder}\t{exp}\t{info}{idfinal}\n')

                #calculated CF and TT
            for file in files:
                if 'mf.tsv' in file:
                    namef = file.split('.')[0]
                    addfile = f'./GEO/{folder}/{file}'
                    analysis(addfile)
                    print (f'READY')
    input("""
                    Expression values have been calculated!!!
                    Check the information in "./ALL/Info_sampes.txt" file, 
                    or the values in each subfolder folder to validadte, 
                    it corresponds to *.xlsx files
          
                    press ENTER to continue""")

def valuesExp():
    folders = os.listdir('./GEO')
    for folder in folders:
        if '.' not in folder:
            files = os.listdir(f'./GEO/{folder}')
            for file in files:
                if 'matrix' in file:
                    outN = open(f'./ALL/{file}', 'w')#to make a copy without info in ALL folder
                    namef = file.split('_')[0]
                    out = open(f'./GEO/{folder}/{namef}.mf.tsv', 'w')
                    file = open(f'./GEO/{folder}/{file}')
                    rfile = file.readline()
                    keyexp = ['!Series_sample_id', '!Sample_title', '\"ID_REF\"', ]
                    keyexpv = ['NONE', 'NONE', 'NONE']
                    infosamples = ''
                    valuesamples = ''
                    while rfile!='':
                        #print (namef, rfile.strip())
                        for k in range(len(keyexp)):
                            if keyexp[k] in rfile:
                                keyexpv[k]='YES'
                                infosamples+=rfile
                    
                        if '\"ID_REF\"' in rfile:
                            outN.write(rfile)#to make a copy without info in ALL folder
                            rfile = file.readline()
                            while '!series_matrix_table_end' not in rfile:
                                outN.write(rfile)#to make a copy without info in ALL folder
                                valuesamples+=rfile
                                rfile = file.readline()
                        rfile = file.readline()
                    
                    print (f"""
                    Quality control for {folder}/{namef}: 
                        include sample ID [{keyexpv[0]}]
                        include sample title [{keyexpv[1]}]
                        include ID_REF [{keyexpv[2]}]
                        """)
                    input ("""
                    If quality control is OK, please define the samples
                    and the kind of sample to each experiment. They will be
                    used to calculate the expression profiles and the
                    comparison among experiments.
                            
                    To do that, press ENTER to continue
                    """)

                    namesam = infosamples.split('\n')[1]
                    classuser = '!UserClassification'
                    for ns in namesam.split('\t')[1:]:
                        cu = input(f'check ({namef}) and classify ({ns}) as control/sample/none - C/S/N: ')
                        classuser+=f'\t{cu}'
                    input ('\nREADY, press ENTER to the next experiment')
                    out.write(f'{infosamples}{classuser}\n{valuesamples}')

    input('\nExperiments have been loaded, Check in each subforlder the file *.mf.tsv to validate, \n\npress ENTER to continue')

def menu():
    os.system('clear')
    print ("""
                    GEO TARGET EXPRESSION PROFILE
                             (GEOTEP)
                            Created by 
                        Fabian Tobar-Tosse
                    ftobar@javerianacali.edu.co
                Pontifica Universidad Javeriana Cali

           
    This program extracts expression profiles of specific genes 
    from expression studies in the GEO database. Its main objective 
    is to identify genes that are co-expressed with target genes in 
    expression experiments, even if they are not the primary focus 
    of the original experiments..

    The next steps enable you to integrate multiple experiments and 
    define specific comparisons based on the samples. Follow the 
    next steps in order.
           
    a. Get data from GEO experiments: 
           
        In this step, you should classify the GEO experiments in folders 
        that define the type of experiments, e.g., folders: ROOTS, STEMS, 
        LEAVES, etc. Besides, you should define controls (C), samples (S), 
        or NONE (N) for each sample in the experiment. This means you can 
        compare only some samples from each experiment, not the same ones 
        explored in the original experiment. Move the expression data into
        the ./GEO folder. Specifically, you need to download the files with
        the extension "_series_matrix.txt" located in the "download data"
        section of each experiment found in GEO database
              
    b. Get Log2 Fold Change and T-test value: 
           
        This step calculates the Log2FC and T-test values for the selected 
        samples (C/S) that the user defined in the last step. In addition, 
        an "./ALL/Info_samples.txt" file is created to save the configuration 
        defined by the user for data comparison.
           
    c. Get genes with similar profiles to the target genes. 
           
        IMPORTANT: For this step, you need a file called “./GEO/genename.txt” 
        that contains the geneID and the gene_name. It is important to define 
        the profiles under a standard name, that could be useful to analyze 
        the profiles. 
           
        You can get this data from the GEO database, mainly the “soft” file
           
    d. Make a comparative matrix based similar profiles:
           
        This step allows to identify common genes among the different experiments,
        and with common profiles. Those experiments are comparable based on the FC
        or T-test scores. 
           
    x. EXIT
           """)
    
    usrO = input('write an option and press ENTER: ')
    if usrO == 'a':
        os.system('clear')
        print ('RUNNING STEP A\n')
        #This folder is necessary to the next step
        folders = os.listdir('./')
        if 'ALL' not in folders:
            os.system('mkdir ALL')
        if 'GEO' not in folders:
            os.system('mkdir GEO')
        #run the funtion
        check1 = input('Do you have the experiments arranged in subfolders into the ./GEO folder? y/n: ')
        if check1 == 'y':
            valuesExp()
            menu()
        else:
            input("""
                  Please make subfolders in the ./GEO folder, 
                  move the files "_series_matrix.txt" into 
                  each one and run again, press ENTER to go back. 
                  """)
            menu()
    elif usrO == 'b':
        os.system('clear')
        print ('RUNNING STEP B\n')
        loadmatrix()
        menu()

    elif usrO == 'c':
        os.system('clear')
        print ('RUNNING STEP C\n')
        listag = input("""
        write the gene(s) id that you like to retrive the expression profiles.
        ex, 262455_at,246312_at,262234_at,260145_at,267136_at,267375_at
        : """)
        listag = listag.split(',')
        fdr = os.listdir('./')
        if 'temp' not in fdr: 
            os.system('mkdir temp')#temporal
        if 'ExpPro' not in fdr: 
            os.system('mkdir ExpPro')#result
        for lg in listag:
            print (f"""
            Identifing the genes with the same expression profile to the ({lg}) gene ...
            """)
            getdata(lg)
            makematrix(lg)
        os.system('rm -rf temp')
        input ("""
        ExpPro folder have been created. It contains the expression
        profiles to each gene in each experiment or comparison, the
        gene_name was included to a better interpretation of the 
        profiles, besides the kind of comparison was written
        with experiment_id. press ENTER to continue
        """)
        menu()
    elif usrO == 'd':
        os.system('clear')
        print ('RUNNING STEP D\n')
        ComparativeMatrix()
        input ("""
        Ready, check the the Comparative_Matrix files in the
        "ExpPro" folder.
        
                    THANKS FOR USING GEOTEP
               
        Press ENTER to exit""")

def warn():
    warnings.filterwarnings("ignore", category=RuntimeWarning)
import os
import pandas as pd
import numpy as np
from scipy import stats
import warnings
warn()
menu()