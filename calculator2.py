import math
from StringIO import StringIO
import numpy
import csv
import operator
import random
from scipy.stats import ks_2samp
from scipy.stats import stats
import itertools
import gzip
import collections
import re

# CONSTANTS would go here, if any - eg, CONSTANT={'value1': num, 'value2':num2, etc.}


def calculate_similarity(variables): 
    """

    This calculates similarity between individual drugs, and groups of other drugs

    Variables is a dictionary passed from server.py
 
    """

# take csv list passed of human genes
    #fixcomplist=re.split(',',variables['Druglist'].splitlines())
    complist=[x.strip() for x in variables['Druglist'].replace('\n',',').split(',')]  
    #complist=[x.strip() for x in fixcomplist] 
    
# read in pairs of human id to worm
    fname="static/orthology_rp.txt"
    #fname="static/toy_orthology_rp.txt"
# read in property dictionary - new format
    matcher_dict={}
    fullgenelist=[]   
    matchedlist=[]
    goodgenes=[]
    
    genepairs=csv.reader(open(fname),delimiter='\t')
    for row in genepairs:
        if row[0]<>'Ensembl Gene ID' : matcher_dict[row[3]]=row[0]
    for comp in complist: 
        if matcher_dict.has_key(comp):
            matchedlist.append(matcher_dict[comp]) 
            goodgenes.append(comp)
    #print (matcher_dict)
    
    fname="static/mart_export.txt.gz"

# read in worm phenotype dictionary - new format
    property_dict={}
    fulldruglist=[]   
    phenolist=[]
    with gzip.open(fname) as gzfile:
        drug_properties=csv.reader(gzfile,delimiter='\t')
        for row in drug_properties:
            templist=row[2].split(' | ')
            tempfulllist=row[3].split(' | ')
            templist2=[]
            tempfulllist2=[]
            for items in templist: 
                if not(items==''): templist2.append(items)
                #### note big problem: in the shorter list, no distinction between observed and NOT!!!
            for items in tempfulllist:
                if not(items=='') and not("Not Observed" in items): tempfulllist2.append(items[10:-47])
            #liststring="|".join(templist)
            property_dict[row[5]]=[row[1],templist2,tempfulllist2,row[4]]
            phenolist.append(row[2].split(' | '))
        for key in property_dict: 
            fulldruglist.append(key) 
    #print(property_dict)
# now output phenotypes
    pheno_list=[]
    fullpheno_list=[]
    goodwormgene=[]
    klist=[]
    fullklist=[]
    for wormgene in matchedlist :
        if property_dict.has_key(wormgene) :
            temppheno=property_dict[wormgene][1]
            tempfullpheno=property_dict[wormgene][2]
          
            #print(templistfixed)
            if temppheno<>[]:
                pheno_list=pheno_list+temppheno
                klist=klist+temppheno
                fullpheno_list=fullpheno_list+tempfullpheno
                fullklist=fullklist+tempfullpheno
            goodwormgene.append(wormgene)
    appearance_counterlist=collections.Counter(klist)
    appearancefull_counterlist=collections.Counter(fullklist)
    #print (klist,appearance_counterlist.most_common())

# now make a list of every phenotype for every gene
    matchedlist2=[]    
    for key in matcher_dict: 
        if matcher_dict.has_key(key):
            matchedlist2.append(matcher_dict[key]) 

# now collate phenotypes
    pheno2_list=[]
    fullpheno2_list=[]
    goodwormgene2=[]
    klist2=[]
    fullklist2=[]
    for wormgene in matchedlist2 :
        if property_dict.has_key(wormgene) :
            temppheno=property_dict[wormgene][1]
            tempfullpheno=property_dict[wormgene][2]
            if not(temppheno==[]):
                pheno2_list=pheno2_list+temppheno
                klist2=klist2+temppheno
                fullpheno2_list=fullpheno2_list+tempfullpheno
                fullklist2=fullklist2+tempfullpheno
            goodwormgene2.append(wormgene)
    appearance_counterlist2=collections.Counter(klist2)
    appearancefull_counterlist2=collections.Counter(fullklist2)
    
    # make a dictionary with pathways
    
    pathway_dict={}
    commonlist=appearancefull_counterlist2.most_common()
    for items in commonlist :
        #pathway_dict[items][1]="nose"
        pathway_dict[items[0]]=[]
    
    # loop through all genes and add them to pathways

    for keys in pathway_dict :
        for genes in matcher_dict :
            #print("key:",keys,"genes",genes)
            #if property_dict.has_key(matcher_dict[genes][2]):
            #print("true",matcher_dict[genes])    
            #print("item:",property_dict[matcher_dict[genes]][2])
            if keys in property_dict[matcher_dict[genes]][2]: 
                pathway_dict[keys].append(genes)
    print (pathway_dict)
    
    return {
        'worm_phenos': pheno_list,
        'goodwormgene':goodwormgene,         
        'goodgenes':goodgenes,
        'pathways_matched':klist,
        'pathway_counts':appearance_counterlist.most_common(),
        'full_pathways':fullklist,
        'full_pathways_counts':appearancefull_counterlist2.most_common()
    }


