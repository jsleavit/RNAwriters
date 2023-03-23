#!/usr/bin/env python3

import pandas as pd
import sys
import os

''' Use Gtrnadb to download NCBI genome data '''

class COLLECT_GENOME_INFO():
    def __init__(self, genomeInfoFile, downloadDirectory):
        self.genomeInfo = pd.read_csv(genomeInfoFile, sep='\t')
        self.downloadDirectory = downloadDirectory

    def downloadNCBIgenomeInfo(self, accs, fileName, downloadDir, downloadOptions = 'genome,gtf,protein,cds,rna'):
        """Use ncbi datasets tool to retrieve specified files"""
            
        if downloadDir[-1] != '/':
            downloadDir += '/'
        
        try:
            os.mkdir(downloadDir)
        except FileExistsError:
            pass
        

        if os.path.exists(fileName):
            pass
        else:
            options = '--include ' + downloadOptions
            os.system(f'datasets download genome accession --inputfile {accs} --no-progressbar {options} --filename {fileName}')

        os.system(f'unzip -d {downloadDir} {fileName}')
        
    
    def getAccessions(self):
        ''' Get list of accessions from genomeInfo file '''

        fileName = f'gtrnadb_genomes.zip'
        tmpAccess = f'genome_accessions_tmp.txt'
        if os.path.exists(fileName):
            pass
        else:
            accessions = self.genomeInfo['Gb_asm_acc'].tolist()
            with open(tmpAccess, 'w') as f:
                for acc in accessions:
                    f.write(acc.strip() + '\n')

        #self.downloadNCBIgenomeInfo(tmpAccess, fileName, self.downloadDirectory)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python3 collect_genome_info.py genomeInfoFile downloadDirectory')
        sys.exit(1)
    COLLECT_GENOME_INFO(sys.argv[1], sys.argv[2]).getAccessions()