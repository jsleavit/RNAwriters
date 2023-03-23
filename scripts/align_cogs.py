#!/usr/bin/env python3

''' align COG sequences '''

import os
import glob
from Bio import AlignIO
class ALIGN_COGS():
    def __init__(self, clusteredFastasList = [],
                    outPath = '/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/aligned_fastas/',
                    threads = 20,
                    overwrite = False):
        
        self.clustered_fastas = clusteredFastasList
        self.out_path = outPath
        self.overwrite = overwrite
        self.threads = threads

    def fasta2sto(self, fastaFile):
        ''' convert fasta file to stockholm format '''
        cog = fastaFile.split('/')[-1].split('.')[0]
        alignment = AlignIO.parse(open(fastaFile), 'fasta')
        AlignIO.write(alignment, fastaFile.replace('.fasta', '.sto'), 'stockholm')
        print(f'alignment for {cog} written to {fastaFile}')

    def alignCogSequences(self): 
        ''' align cog sequences '''
        for path in self.clustered_fastas:
            cog = path.split('/')[-1].split('.')[0]
            out_file = self.out_path + f'{cog}.aligned.clustered.fasta'
            if os.path.exists(out_file) and not self.overwrite:
                print(f'alignment file {out_file} exists')
                pass
            else:
                os.system(f'mafft --thread {self.threads} --auto {path} > {out_file}')

                if not os.path.exists(out_file.replace('.fasta', '.sto')):
                    
                    self.fasta2sto(out_file)
        
    
if __name__ == '__main__':

    # align COGs
    clusteredFastasList = [x for x in glob.glob('/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/clustered_fastas/*.clustered.fasta') if 'clustered' in x]
    ALIGN_COGS(clusteredFastasList).alignCogSequences()
