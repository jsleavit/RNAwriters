#!/bin/usr/env python3

import os, sys
import glob
import pandas as pd

## HMM PARSE STATS
############################################################################################################
import pandas as pd
class HMM_PARSE():
    def __init__(self, hmmsearchResultsList = [],
                    outPath = '/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/',
                    overwrite = False):
        
        self.hmmsearch_results = hmmsearchResultsList
        self.out_path_root = outPath
        self.overwrite = overwrite
        self.genomeInfo = pd.read_csv('/projects/lowelab/users/jsleavit/git_repos/data/Gtrnadb_prokaryotic_info_with_full_taxonomy_updated_asm_accessions.tsv', sep = '\t', index_col = 0)
        
    
    def makeDirectory(self, directory):
        ''' make a directory if it does not exist '''
        if not os.path.exists(directory):
            os.makedirs(directory)
        if os.path.exists(directory):
            return True

    def buildResultPath(self):
        ''' build path to hmmsearch results '''

        # write hmmsearch_result paths to file for later use
        self.makeDirectory(self.out_path_root + 'hmmsearch_stats/')
        if os.path.exists(self.out_path_root + 'hmmsearch_stats/hmmsearch_paths.txt') and not self.overwrite:
            print(f'assembly list {self.out_path_root + "hmmsearch_stats/hmmsearch_paths.txt"} exists')
            pass
        else:
            with open(self.out_path_root + 'hmmsearch_stats/hmmsearch_paths.txt', 'w') as f:
                for hmmsearch_result in self.hmmsearch_results:
                    f.write(hmmsearch_result + '\n')

    
    def checkAssemblies(self):
        ''' check what assemblies are in the hmmsearch results '''

        if os.path.exists(self.out_path_root + 'hmmsearch_stats/hmmsearch_paths.txt') and not self.overwrite:
            pass
        else:
            self.buildResultPath()

        # check assemblies that intersect with gtrnadb
        assemblies = []
        for hmmsearch_result in open(self.out_path_root + 'hmmsearch_stats/hmmsearch_paths.txt', 'r'):
            assemblies.append(hmmsearch_result.split('/')[-1].split('.hmmsearch.out')[0])
        
        domain_counts = {}
        for genome_acc in self.genomeInfo['Gb_asm_acc']:
            if genome_acc in assemblies:
                domain = self.genomeInfo.loc[self.genomeInfo['Gb_asm_acc'] == genome_acc, 'Domain'].values[0]
                domain_counts[domain] = domain_counts.get(domain, 0) + 1
        
        gtRNAdb_counts = self.genomeInfo['Domain'].value_counts()
        
        print([f'{domain}: {domain_counts[domain]} ({round(domain_counts[domain]/gtRNAdb_counts[domain]*100, 2)}%)' for domain in domain_counts.keys() if domain in gtRNAdb_counts.keys()])
        
    def parseHMMsearch(self):
        ''' parse hmmsearch results '''

        if os.path.exists(self.out_path_root + 'hmmsearch_stats/hmmsearch_paths.txt') and not self.overwrite:
            pass
        else:
            self.buildResultPath()

        # parse hmmsearch results
    
        for hmmsearch_result in open(self.out_path_root + 'hmmsearch_stats/hmmsearch_paths.txt', 'r'):
            genome = hmmsearch_result.split('/')[-1].split('.hmmsearch.out')[0]
            phylum = self.genomeInfo.loc[self.genomeInfo['Gb_asm_acc'] == genome, 'phylum'].values[0]
            family = self.genomeInfo.loc[self.genomeInfo['Gb_asm_acc'] == genome, 'family'].values[0]
            species = self.genomeInfo.loc[self.genomeInfo['Gb_asm_acc'] == genome, 'species'].values[0]

            if self.makeDirectory(self.out_path_root + f'hmmsearch_stats/{genome}/'):
                pass
            else:
                self.makeDirectory(self.out_path_root + f'hmmsearch_stats/{genome}/')
            
            cog_genes = {}
            with open(hmmsearch_result.strip(), 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    if 'clustered' in line.split()[2]:
                        cog = line.split()[2].split('.')[0]
                        gene = line.split()[0]
                        e_val = line.split()[4]
                        product = line.split()[-1]
                        cog_genes[cog] = cog_genes.get(cog, []) + [(gene, e_val, product)]
                
            
            for cog in cog_genes.keys():
                cog_gene_count = len(cog_genes[cog])
                out_file = self.out_path_root + f'hmmsearch_stats/{genome}/{genome}.{str(cog_gene_count)}.{cog}.hmmsearch.out'
                if os.path.exists(out_file) and not self.overwrite:
                    print(f'hmmsearch file {out_file} exists')
                    pass
                else:
                    with open(out_file, 'w') as f:
                        for gene, e_val, product in cog_genes[cog]:
                            to_write = [str(x).strip() for x in [genome, cog, cog_gene_count, gene, e_val, product, phylum, family, species]]
                            f.write('\t'.join(to_write) + '\n')
                            #f.write(f'{genome}\t{cog}\t{gene}\t{e_val}\t{product}\t{phylum}\t{family}\t{species}\n')
                            #f.write(f'{family}{species}\t{cog}{gene}\t{e_val}\t{product}\n')
if __name__ == '__main__':

    # move hmmsearch results
    #hmmsearch_results = glob.glob('/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/hmmsearch_genomes/*.out')
    #HMM_STATS(hmmsearchResultsList = hmmsearch_results).checkAssemblies()
    #HMM_PARSE().checkAssemblies()
    HMM_PARSE().parseHMMsearch()