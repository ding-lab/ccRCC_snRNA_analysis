#jupyter notebook
#conda install pandas
import argparse
import sys
import pandas as pd
import numpy as np
import scipy
from scipy.stats import chisquare
import warnings
import statsmodels.api

def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("-p", "--protein_exp", help="Your protein expression matrix downloaded from HPA")
    parser.add_argument("-o", "--output", help="Your destination output file.")
    parser.add_argument("-t", "--tissue_type", help="Tissue type summary metadata")
    parser.add_argument("-g", "--gene_list",help="Data frame of gene list")
    options = parser.parse_args(args)
    return options

options = getOptions(sys.argv[1:])


#pro_df=pd.read_csv("/diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Data/HPA_normal_tissue.tsv",delimiter='\t')
pro_df=pd.read_csv(options.protein_exp,delimiter='\t')
pro_df.rename(columns={'Gene name': 'Gene_name',"Cell type":"Cell_type"}, inplace=True)
pro_df=pro_df.replace(' ', '_', regex=True)

#self-made reference table to convert detailed tissue types in to general tissue types shown on the webpage
#tissue_type_df=pd.read_csv("/diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Data/HPA_Tissue_type_matching.txt",delimiter='\t')
tissue_type_df=pd.read_csv(options.tissue_type,delimiter='\t')
#out_path="/diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Results/V8/Union_discovery"
out_path=options.output

#function 
warnings.simplefilter("ignore")

#function to 1) parse the protein expression table from PROTEIN ATLAS 2) test if a certain gene is significantly expressed in certain tissues
def pro_exp_subset(gene):
    pro_df_subset=pro_df.loc[pro_df['Gene_name']==gene]
    pro_df_subset['general_tissue']=pro_df_subset['Tissue'].map(tissue_type_df.set_index('tissue')['general_tissue'])

    pro_df_subset['Level'] = np.where((pro_df_subset.Level == 'Low'),1,pro_df_subset.Level)
    pro_df_subset['Level'] = np.where((pro_df_subset.Level == 'Medium'),2,pro_df_subset.Level)
    pro_df_subset['Level'] = np.where((pro_df_subset.Level == 'High'),3,pro_df_subset.Level)
    pro_df_subset['Level'] = np.where((pro_df_subset.Level == 'Not_detected'),0,pro_df_subset.Level)

    tmp = pro_df_subset.groupby('general_tissue')['Level'].max().reset_index()
    tmp['Gene']=pro_df_subset.iloc[0,0]
    tmp['Gene_name']=gene
    tmp=tmp.reindex(columns=['Gene','Gene_name','Level','general_tissue'])
    
    #calculate chi-square pvalue. Assuming the distribution if uniform distribution
    colnames=[0,1,2,3]
    freq_table=pd.DataFrame(index=["data","obs_expected"],columns=colnames)
    freq_table.loc[['obs_expected']]=tmp.shape[0]/4
    obs=[]
    for col in colnames:
        freq_table.loc[["data"],[col]]=tmp[tmp['Level'] == col].shape[0]
        obs.append(tmp[tmp['Level'] == col].shape[0])
    pval=chisquare(obs)[1]
    
    max_exp=tmp['Level'].max()
    tmp_high=tmp.loc[tmp['Level']==max_exp]

    return(gene,pval,tmp_high);


#gene_list=pd.read_csv("/diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Results/V8/Union_discovery/Gene_list_of_DE_genes_filtered_surface_3DB.txt",header=None,delimiter='\t')
gene_list=pd.read_csv(options.gene_list,header=None,delimiter='\t')
genes=gene_list[0].to_list()
len(genes)

pro_exp_df=pd.DataFrame(columns=["Gene","Gene_name","Level","general_tissue"])
pval_table=pd.DataFrame(index=genes,columns=["gene","pval"])
for gene in genes:
    if gene in pro_df['Gene_name'].tolist():
        pval_table.loc[[gene],"gene"]=pro_exp_subset(gene)[0]
        pval_table.loc[[gene],"pval"]=pro_exp_subset(gene)[1]  
        pro_exp_df=pd.concat([pro_exp_df,pro_exp_subset(gene)[2]],ignore_index=True)
    else:
        pval_table.loc[[gene],"gene"]="NaN"
        pval_table.loc[[gene],"pval"]="NaN"
        
pval_table=pval_table.drop(pval_table[pval_table['pval'] == "NaN" ].index) #remove genes that are not found in pro database

#split table based on whether pvalue is NaN or not, NAN means no protein data in the HPA database
pval_table_ori=pval_table
row_has_NaN=pval_table_ori['pval'].isnull()
pval_table_with_NaN = pval_table_ori[row_has_NaN]
pval_table_without_NaN = pval_table_ori[~pval_table_ori.index.isin(pval_table_with_NaN.index)]

################################################
### FDR fdr_bh
#################################################
#conda install -c conda-forge statsmodels
#from statsmodels.stats._knockoff import RegressionFDR
pval_table=pval_table_without_NaN
pvals=pval_table['pval'].tolist()
pval_table['FDR']=statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.05, method='i', is_sorted=False)[1]
pval_table['FDR_filter']=statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.05, method='i', is_sorted=False)[0]

#add 2 'FDR' columns to keep the same format as pval_table_without_NaN and then combine the rows with/without pvalues
pval_table_with_NaN['FDR']="NaN"
pval_table_with_NaN['FDR_filter']="NaN"
pval_table_comb=pval_table.append(pval_table_with_NaN)
pval_table=pval_table_comb

#only keep genes if they pass FDR_fitler 
keep_genes=pval_table[pval_table['FDR_filter']==True]['gene'].to_list()
pro_exp_df=pro_exp_df[pro_exp_df['Gene_name'].isin(keep_genes)]

#rename colnames
pro_exp_df= pro_exp_df.rename(columns = {'Gene':'GENE','Gene_name':'GENE_SYMBOL','Level':'EXP_LEVEL','general_tissue':'GENERAL_TISSUE_TYPE'})
pval_table=pval_table.rename(columns={'gene':'GENE','pval':'P-VALUE','FDR_filter':'FDR_FILTER'})

#save df, significant proteins 
pro_exp_df.to_csv(out_path+"/HPA_Protein_expression_for_tissue_specific_DEG.txt",index=False,sep="\t")
#pval_table contains all the input genes even the genes that could not be found in HPA database, which is labeled as NA in the end of the this file.
pval_table.to_csv(out_path+"/HPA_Protein_significance_test_pval_of_all_input_genes.txt",index=False,sep="\t")
