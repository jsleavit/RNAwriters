#!/usr/bin/env python3

import os
import sys
import glob

# HMM_COGS
############################################################################################################
class HMM_COGS():
    def __init__(self, alignedFastasList = [], genomeDir = '/projects/lowelab/users/jsleavit/git_repos/data/ncbi/ncbi_dataset/data/',
                    outPath = '/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/',
                    threads = 25,
                    overwrite = False):
        
        self.aligned_fastas = alignedFastasList
        self.genome_dir = genomeDir
        self.out_path_root = outPath
        self.overwrite = overwrite
        self.threads = threads
        self.cog_hmms = []

    def buildHMM(self):
        ''' build hmm profile from stockholm alignment '''

        self.out_path_root = self.out_path_root + 'hmmbuilds/'
        if not os.path.exists(self.out_path_root):
            os.mkdir(self.out_path_root)

        for aln in self.aligned_fastas:
            cog = aln.split('/')[-1].split('.')[0]
            hmm_file = self.out_path + 'hmmbuilds/' + f'{cog}.hmm'
            if os.path.exists(hmm_file) and not self.overwrite:
                print(f'hmm file {hmm_file} exists')
                pass
            else:
                print(f'building hmm for {cog}')
                os.system(f'hmmbuild --cpu {self.threads} {hmm_file} {aln} > {hmm_file.replace(".hmm", ".log")}')
                #print(f'hmm for {cog} written to {hmm_file}')

    def buildDatabase(self):
        ''' concatenate all hmm profiles into a single database '''

        hmm_files = glob.glob([x for x in self.out_path + 'hmmbuilds/*.hmm' if 'db' not in x])
        hmm_db = self.out_path + 'hmmbuilds/cog_hmm_db.hmm'
        if os.path.exists(hmm_db) and not self.overwrite:
            print(f'hmm database {hmm_db} exists')
            pass
        else:
            print(f'building hmm database')
            os.system(f'cat {" ".join(hmm_files)} > {hmm_db}')
    
        # press the database
        os.system(f'hmmpress {hmm_db}')


    def searchHMM(self):
        ''' search hmm profile against genomes '''
        
        hmm_db = self.out_path + 'hmmbuilds/cog_hmm_db.hmm'

        assembly_proteins = []
        for assembly in os.listdir(self.genome_dir):
            if os.path.isdir(self.genome_dir+assembly):
                for file in os.listdir(self.genome_dir+assembly):
                    if file == 'protein.faa':
                        assembly_proteins.append(self.genome_dir+assembly+'/'+file)

        # make a new directory for the hmmsearch results
        self.out_path = self.out_path + 'hmmsearch_genomes/'
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)
        
        # search hmm database against all genomes
        for assembly in assembly_proteins:
            genome = assembly.split('/')[-2]
            out_file = self.out_path + f'{genome}.hmmsearch.out'
            if os.path.exists(out_file) and not self.overwrite:
                print(f'hmmsearch file {out_file} exists')
                pass
            else:
                os.system(f'hmmsearch --tblout {out_file} --cpu {self.threads} {hmm_db} {assembly} > {out_file.replace(".out", ".log")}')
                
if __name__ == '__main__':
    # aligned COGs
    #alignedFastasList = [x for x in glob.glob('/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/aligned_fastas/*.sto') if 'sto' in x]
    #HMM_COGS(alignedFastasList = alignedFastasList).buildHMM()

    # build hmm database
    #HMM_COGS().buildDatabase()

    # search hmm database against genomes
    HMM_COGS().searchHMM()