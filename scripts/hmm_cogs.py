#!/usr/bin/env python3

import os
import sys
import glob

# HMM_COGS
############################################################################################################
class HMM_COGS():
    def __init__(self, alignedFastasList = [], genomeDir = '',
                    outPath = '/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/hmms/',
                    threads = 8,
                    overwrite = False):
        
        self.aligned_fastas = alignedFastasList
        self.genome_dir = genomeDir
        self.out_path = outPath
        self.overwrite = overwrite
        self.threads = threads
        self.cog_hmms = []

    def buildHMM(self):
        ''' build hmm profile from stockholm alignment '''

        for aln in self.aligned_fastas:
            cog = aln.split('/')[-1].split('.')[0]
            hmm_file = self.out_path + f'{cog}.hmm'
            if os.path.exists(hmm_file) and not self.overwrite:
                print(f'hmm file {hmm_file} exists')
                pass
            else:
                os.system(f'hmmbuild --cpu {self.threads} {hmm_file} {aln} > {hmm_file.replace(".hmm", ".log")}')
                print(f'hmm for {cog} written to {hmm_file}')

    def searchHMM(self):
        ''' search hmm profile against genomes '''
        
        
        cog_hmms = glob.glob(self.out_path + '*.hmm')

        assembly_proteins = []
        for assembly in os.listdir(self.genome_dir):
            if os.path.isdir(self.genome_dir+assembly):
                for file in os.listdir(self.genome_dir+assembly):
                    if file == 'protein.faa':
                        assembly_proteins.append(self.genome_dir+assembly+'/'+file)
            
        for hmm in cog_hmms:
            cog = hmm.split('/')[-1].split('.')[0]
            for protein in assembly_proteins:
                genome = protein.split('/')[-2]
                out_file = self.out_path + f'{cog}.{genome}.hmmsearch.out'
                if os.path.exists(out_file) and not self.overwrite:
                    print(f'hmmsearch file {out_file} exists')
                    pass
                else:
                    os.system(f'hmmsearch --cpu {self.threads} {hmm} {protein} > {out_file}')
                    print(f'hmmsearch for {cog} written to {out_file}')


    
if __name__ == '__main__':
        # align COGs
    alignedFastasList = [x for x in glob.glob('/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/aligned_fastas/*.sto') if 'sto' in x]
    #HMM_COGS(alignedFastasList = alignedFastasList).buildHMM()
    HMM_COGS(alignedFastasList = alignedFastasList, genomeDir = '/projects/lowelab/users/jsleavit/git_repos/data/ncbi/ncbi_dataset/data/').searchHMM()