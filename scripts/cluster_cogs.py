#!/bin/usr/env python3

import os, sys
import glob
import pandas as pd
import io

# CLUSTER_COGS
############################################################################################################
import gzip
import io
import numpy as np
class CLUSTER_COGS():
    def __init__(self, cogFastaPathList = [],
                  tempPath = '/projects/lowelab/users/jsleavit/git_repos/temp/',
                    outPath = '/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/clustered_fastas/',
                    clusterIdentity = 0.8,
                    randomSelection = 500,
                    coverageMode = 0,
                    minCoverage = 0.8,
                    overwrite = False,
                    subset = True):
        
        self.cog_fasta_path = cogFastaPathList
        self.temp_path = tempPath
        self.out_path = outPath
        self.cluster_identity = clusterIdentity
        self.random_selection = randomSelection
        self.coverage_mode = coverageMode
        self.min_coverage = minCoverage
        self.overwrite = overwrite
        self.subset = subset


    def readFasta(self, inFile = ''):
        ''' read in fasta file and return dictionary of sequences '''

        seqs = {}
        if inFile.split('.')[-1] == 'gz':
            with gzip.open(inFile, 'rb') as inF:
                with io.TextIOWrapper(inF, encoding='utf-8') as decoder:
                    for seq in decoder.read().split('>')[1:]:
                        seqs[seq.split('\n')[0]] = ''.join(seq.split('\n')[1:])
            inF.close()
        else:
            with open(inFile, 'r') as inF:
                for seq in inF.read().split('>')[1:]:

                    seqs[seq.split('\n')[0]] = ''.join(seq.split('\n')[1:])

                inF.close()
        return seqs
    
    def writeFasta(self, data, fName, blockSize = 60):
   
        ''' write dictionary of sequences to fasta file '''
        with open(fName, 'w') as f:
            for seq in data.keys():
                f.write('>{0}\n'.format(seq))
                s = data[seq]
                for x in range(0, len(s), blockSize):
                    f.write(s[x:x+60] + '\n')
            f.close()

    
    def cluster(self, fastas, cog):
        ''' use mmseqs2 to cluster sequences '''

        fasta_string = ' '.join(fastas)

        output_path = f'{self.out_path}{cog}'

        cluster_file = f'{self.out_path}{cog}_cluster.tsv'

        if os.path.exists(cluster_file) and not self.overwrite:
            print(f'cluster file {cluster_file} exists')
            pass
        else:
            print(f'clustering {cog}')
            os.system(f'mmseqs easy-cluster {fasta_string} {output_path} {self.temp_path} --min-seq-id {self.cluster_identity} --cov-mode {self.coverage_mode} -c {self.min_coverage}')

        return cluster_file
    
    def subsetRepresentativeSequences(self, clusterFile, cog, cogSeqs):
        ''' subset representative sequences from cluster file '''

        rep_seq_fasta = self.out_path + f'{cog}.clustered.fasta'
        rep_seqs = pd.read_csv(clusterFile, sep = '\t', header = None, index_col = 1)
        print(len(set(rep_seqs.index.values))),'reduced to', len(set(rep_seqs[0]))

        rep_dict = {}
        if len(set(rep_seqs[0])) > self.random_selection:
            rep_seqs = rep_seqs.sample(self.random_selection)
            subset = np.array(list(set(rep_seqs[0])))
            rep_dict = {x: cogSeqs[x] for x in subset[:self.random_selection]}
        else:
            rep_dict = {x: cogSeqs[x] for x in set(rep_seqs[0])}

        print(f'writing {len(rep_dict)} representative sequences to {rep_seq_fasta}')
        self.writeFasta(rep_dict, rep_seq_fasta)
    
    def clusterCogSequences(self):
        ''' get all sequences from COGs and cluster them '''

        for path in self.cog_fasta_path:
            cog = path.split('/')[-1].split('.')[0]
            if os.path.exists(path.replace('.gz', '')):
                cog_seqs = self.readFasta(inFile = path)
                cog_seqs = {x.split(' ')[0]:cog_seqs[x] for x in cog_seqs}
                cluster_file = self.cluster([path], cog = cog)

                if self.subset:
                    self.subsetRepresentativeSequences(cluster_file, cog, cog_seqs)

if __name__ == '__main__':
    # download COGs
    out_path = '/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/fastas/'
    #cogs_list = [x.strip('\n').split('\t')[0] for x in list(open('/projects/lowelab/users/jsleavit/git_repos/data/combined_trna_mod_database.tsv').readlines())][1:]
    #downloadedCogs = GET_COGs().downloadCOGs(out_path, cogs_list)

    # cluster COGs
    #CLUSTER_COGS(cogFastaPathList = downloadedCogs).clusterCogSequences()