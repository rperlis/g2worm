import csv
import gzip
import collections
import os
import random
from scipy.stats.mstats import mquantiles
     
# CONSTANTS would go here, if any - eg, CONSTANT={'value1': num, 'value2':num2, etc.}
FNAME_WORM_PHENOS=os.path.dirname(__file__)+"/static/mart_export.txt.gz" 
FNAME_GENE_MATCHER=os.path.dirname(__file__)+"/static/orthology_rp.txt"

property_dict={}
phenolist=[]
with gzip.open(FNAME_WORM_PHENOS) as gzfile:
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
            if not(items=='') and not("Not Observed" in items): 
                temp_parse=items[10:-47]
                    #print(temp_parse,"end=",temp_parse[-3:])
                if temp_parse.endswith(': ex'): temp_parse=temp_parse[:-4]
                if temp_parse.endswith(': e'): temp_parse=temp_parse[:-3]
                if temp_parse.endswith(': '): temp_parse=temp_parse[:-2]
                if temp_parse.endswith(':'): temp_parse=temp_parse[:-1]
                tempfulllist2.append(temp_parse)
            #liststring="|".join(templist)
        property_dict[row[5]]=[row[1],templist2,tempfulllist2,row[4]]
        phenolist.append(row[2].split(' | '))
            
            
def check_genelist(variables): 
    """

    Check a list of genes against worm phenotypes

    Variables is a dictionary passed from server.py
 
    """

# take csv list passed of human genes
    #fixcomplist=re.split(',',variables['Druglist'].splitlines())
    complist=[x.strip() for x in variables['Druglist'].replace('\n',',').split(',')]  
    #complist=[x.strip() for x in fixcomplist] 
    
# read in pairs of human id to worm

    #fname="static/toy_orthology_rp.txt"
# read in property dictionary - new format
    matcher_dict={}
    fullgenelist=[]   
    matchedlist=[]
    goodgenes=[]
    matchedworm_to_human={}
    
    genepairs=csv.reader(open(FNAME_GENE_MATCHER),delimiter='\t')
    for row in genepairs:
        if row[0]<>'Ensembl Gene ID' : matcher_dict[row[3]]=row[0]
    for comp in complist: 
        if matcher_dict.has_key(comp):
            matchedlist.append(matcher_dict[comp]) 
            goodgenes.append(comp)
            matchedworm_to_human[matcher_dict[comp]]=comp
    print(matchedworm_to_human)
    #print (matcher_dict)
    
#    fname="static/mart_export.txt.gz"

# read in worm phenotype dictionary - new format
# moved to not be a function so gets compiled when server run, just for speediness - thanks Brett!

# now output phenotypes
    pheno_list=[]
    fullpheno_list=[]
    goodwormgene=[]
    klist=[]
    fullklist=[]
    gene_to_pathway={}
    pathway_with_genes={}
    
    for wormgene in matchedlist :
        if property_dict.has_key(wormgene) :
            temppheno=property_dict[wormgene][1]
            tempfullpheno=property_dict[wormgene][2]
            gene_to_pathway[matchedworm_to_human[wormgene]]=tempfullpheno   # assign phenos to human genes

            #print(templistfixed)
            if temppheno<>[]:
                pheno_list=pheno_list+temppheno
                klist=klist+temppheno
                fullpheno_list=fullpheno_list+tempfullpheno
                fullklist=fullklist+tempfullpheno

            #experimental - tabulate number of hits per pathway
                for phenos in tempfullpheno:
                    if pathway_with_genes.has_key(phenos) :
                        if not(matchedworm_to_human[wormgene] in pathway_with_genes[phenos]):
                            pathway_with_genes[phenos].append(matchedworm_to_human[wormgene])
                    else : pathway_with_genes[phenos]=[matchedworm_to_human[wormgene]]
            goodwormgene.append(wormgene)       
    appearance_counterlist=collections.Counter(klist)
    appearancefull_counterlist=collections.Counter(fullklist)
    print ("paths:",pathway_with_genes)

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
    
    
    # return 
    return {
        'worm_phenos': pheno_list,
        'goodwormgene':goodwormgene,         
        'goodgenes':goodgenes,
        'pathways_matched':klist,
        'pathway_counts':appearance_counterlist.most_common(),
        'full_pathways':fullklist,
        'full_pathways_counts':appearancefull_counterlist.most_common(),
        'pathway_with_genes':make_table_list(pathway_with_genes,'pathway','genes'),
        'gene_to_pathway':make_table_list(gene_to_pathway,'gene','pathway')
    }



def make_pathways(variables): 
    """

    This outputs pathway files

    Variables is a dictionary passed from server.py
 
    """

    
# read in pairs of human id to worm
    fname="static/orthology_rp.txt"
    #fname="static/toy_orthology_rp.txt"
# read in property dictionary - new format
    matcher_dict={}
    goodgenes=[]
    
    genepairs=csv.reader(open(fname),delimiter='\t')
    for row in genepairs:
        if row[0]<>'Ensembl Gene ID' : matcher_dict[row[3]]=row[0]
    
    goodgenes=["Not_applicable - creating pathway files"]
    #print (matcher_dict)
    
    fname="static/mart_export.txt.gz"

# read in worm phenotype dictionary - new format
    property_dict={}
 
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
                if not(items=='') and not("Not Observed" in items): 
                    temp_parse=items[10:-47]
                    #print(temp_parse,"end=",temp_parse[-3:])
                    if temp_parse.endswith(': ex'): temp_parse=temp_parse[:-4]
                    if temp_parse.endswith(': e'): temp_parse=temp_parse[:-3]
                    if temp_parse.endswith(': '): temp_parse=temp_parse[:-2]
                    if temp_parse.endswith(':'): temp_parse=temp_parse[:-1]
                    tempfulllist2.append(temp_parse)
            #liststring="|".join(templist)
            property_dict[row[5]]=[row[1],templist2,tempfulllist2,row[4]]
            phenolist.append(row[2].split(' | '))
            #fulldruglist.append(key) 
 
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

    ####################
    # make saved permutation file
    ####################
    pathway_counts={}
    pathway_maincounter={}
    
    for path in pathway_dict:
        pathway_maincounter[path]=[]
    for i in range(5):
        randomgenelist=random.sample(list(matcher_dict.keys()),6)
        print randomgenelist
        for path in pathway_dict:
            pathway_counts[path]=0
            for gene in randomgenelist:
                if gene in pathway_dict[path]:
                    pathway_counts[path]=pathway_counts[path]+1
            pathway_maincounter[path].append(pathway_counts[path])
            print(pathway_counts[path])
   
    pathway_pvals={}
    for path in pathway_dict: 
        pathway_pvals[path]=mquantiles(pathway_dict[path],[0.95,0.99]).tolist()
        print pathway_pvals[path]
        
    ####################
    # write files
    ####################
    inrich_file=file("static/g2worm_inrich","w")
    inrich_file.write("Gene_name\tPath_name\tAnnotation\n")
    for human_genes in matcher_dict :
        for phenotypes in property_dict[matcher_dict[human_genes]][2]:
            inrich_file.write(human_genes)
            inrich_file.write("\t")
            inrich_file.write(phenotypes.replace(" ","_"))
            inrich_file.write("\t")
            inrich_file.write(phenotypes.replace(" ","_"))
            inrich_file.write("\n")
    inrich_file.close()
    
    ######################
    # write set-list hg18 for plink
    ######################
    # first read hg18 data
    # read in pairs of human id to worm
    fname="static/glist-hg18"
    #fname="static/toy_orthology_rp.txt"
# read in property dictionary - new format
    hg18_lookup={}
    gene_dict={}
    
    gene_table=csv.reader(open(fname),delimiter=' ')
    for row in gene_table:
        gene_dict[row[3]]=row[0]+"\t"+row[1]+"\t"+row[2]
# now output list of genes
    plink_file=file("static/g2worm_plinkset","w")
    # plink_file.write("Gene_name\tPath_name\tAnnotation\n") ---> NO HEADER
    for human_genes in matcher_dict :
        for phenotypes in property_dict[matcher_dict[human_genes]][2]:
            if gene_dict.has_key(human_genes):
                plink_file.write(gene_dict[human_genes])
                plink_file.write("\t")
                plink_file.write(human_genes)
                plink_file.write("\t")            
                plink_file.write(phenotypes.replace(" ","_"))
                plink_file.write("\n")
            #else : print("unmatched_hg18_gene:",human_genes)    
    plink_file.close()
    
    # return 
    return {
        'worm_phenos': "Not_applicable",
        'goodwormgene':"Not_applicable",         
        'goodgenes':"Not_applicable",
        'pathways_matched':"Not_applicable",
        'pathway_counts':"Not_applicable",
        'full_pathways':"Not_applicable",
        'full_pathways_counts':"Not_applicable",
        'pathway_with_genes':"Not_applicable",
        'gene_to_pathway':"Not_applicable",
    }

def make_table(any_dictionary,col1,col2,showfreq=1):
    htmltext='<table class="table table-condensed table-bordered table-striped"><thead><th style="width:20%">'+col1+'</th><th>'+col2+'</th></thead>'
    htmltext=htmltext+'<tbody>'
    if any_dictionary:
        for key1 in any_dictionary:
            htmltext=htmltext+'<tr><td>'+key1+'</td><td>'
            for key2 in any_dictionary[key1]:
                if showfreq==1:
                    htmltext=htmltext+key2+' ('+str(any_dictionary[key1][key2]*100)+'%), '
                else:
                    htmltext=htmltext+key2+', '
            htmltext=htmltext[:-2]+'</td></tr>'    
    else: htmltext=htmltext+'<tr><td>None</td><td> </td></tr>'
    htmltext=htmltext+'</tbody></table>'
    return htmltext

def make_table_list(any_dictionary,col1,col2):
    htmltext='<table class="table table-condensed table-bordered table-striped"><thead><th style="width:40%">'+col1+'</th><th>'+col2+'</th></thead>'
    htmltext=htmltext+'<tbody>'
    if any_dictionary:
        for key1 in any_dictionary:
            htmltext=htmltext+'<tr><td>'+key1+'</td><td>'+', '.join(any_dictionary[key1])+'</td></tr>'
    else: htmltext=htmltext+'<tr><td>None</td><td> </td></tr>'
    htmltext=htmltext+'</tbody></table>'
    return htmltext
