#!/usr/bin/env python3

import pandas as pd
import sys
import os
import re

''' Use Gtrnadb to download NCBI genome data '''

class COLLECT_GENOME_INFO():
    def __init__(self, genomeInfoFile, downloadDirectory):
        #self.genomeInfo = pd.read_csv(genomeInfoFile, sep='\t')
        self.genomeInfo = genomeInfoFile
        self.downloadDirectory = downloadDirectory
        self.genomePath = f'{downloadDirectory}ncbi/ncbi_dataset/data/'
        

      #dsets download genome accession 
    def downloadNCBIgenomeInfo(self, accs, fileName, downloadDir, downloadOptions = 'genome,rna,protein,cds,gff3,gtf,gbff,seq-report', source = 'RefSeq'):
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
            options = '--include ' + downloadOptions # genome,gtf,protein,cds,rna 
            os.system(f'datasets download genome accession --inputfile {accs} --no-progressbar {options} --filename {downloadDir}{fileName} --assembly-source {source}')
            os.system(f'unzip {downloadDir}{fileName} -d {downloadDir}ncbi')

    def downloadrRNA(self, genomeRNA_dict):
        """Download rRNA sequences from NCBI using datasets tool"""

        # download rRNA
        for genomePath in genomeRNA_dict:
            if genomeRNA_dict[genomePath]['23S'] != '':
                os.system(f'datasets download gene gene-id {genomeRNA_dict[genomePath]["23S"]} --filename {genomePath}23S.zip --include gene,rna --no-progressbar')
                os.system(f'unzip {genomePath}23S.zip -d {genomePath}23S/')
                os.system(f'mv {genomePath}23S/ncbi_dataset/data/gene.fna {genomePath}23S.fasta')

            if genomeRNA_dict[genomePath]['16S'] != '':
                os.system(f'datasets download gene gene-id {genomeRNA_dict[genomePath]["16S"]} --filename {genomePath}16S.zip --include gene,rna --no-progressbar')
                os.system(f'unzip {genomePath}16S.zip -d {genomePath}16S/')
                os.system(f'mv {genomePath}16S/ncbi_dataset/data/gene.fna {genomePath}16S.fasta')

            if genomeRNA_dict[genomePath]['5S'] != '':
                os.system(f'datasets download gene gene-id {genomeRNA_dict[genomePath]["5S"]} --filename {genomePath}5S.zip --include gene,rna --no-progressbar')
                os.system(f'unzip {genomePath}5S.zip -d {genomePath}5S/')
                os.system(f'mv {genomePath}5S/ncbi_dataset/data/gene.fna {genomePath}5S.fasta')
   

    def getrRNA(self):
        """Read genome annotation file and get the GeneID for 23S ribosomal RNA then write to file and use datasets to download the sequence"""

        genomeRNAID = {}
        with open(self.genomeInfo, 'r') as f:
            accesions = f.readlines()

        
        for acc in accesions:

            # check if genome has been downloaded
            if not os.path.exists(self.genomePath + acc.strip()):
                continue
            
            # check if genomic.gff file exists
            if not os.path.exists(self.genomePath + acc.strip() + '/genomic.gff'):
                continue

            genomePath = self.genomePath + acc.strip() + '/'
            annotationFile = self.genomePath + acc.strip() + '/genomic.gff'

            genomeRNAID[genomePath] = {'16S':'', '23S':'', '5S':''}

            # Get GeneID for 23S rRNA in the genome annotation file
            with open(annotationFile, 'r') as f:
                for line in f:

                    # search for the product name in the line using regular expressions
                    if re.search(r'\b23S ribosomal RNA\b', line):
                        # get the GeneID
                        #tokens = line.split()[8].split(';')

                        # use regular expressions to get the GeneID
                        geneID = re.search(r'GeneID:(\d+)', line)
                        if geneID:
                            value = geneID.group(1)
                            genomeRNAID[genomePath]['23S'] = value
                        

                    if re.search(r'\b16S ribosomal RNA\b', line):
                        # get the GeneID
 
                        # use regular expressions to get the GeneID
                        geneID = re.search(r'GeneID:(\d+)', line)
                        if geneID:
                            value = geneID.group(1)
                            genomeRNAID[genomePath]['16S'] = value
    

                    if re.search(r'\b5S ribosomal RNA\b', line):
                        # get the GeneID
                        #tokens = line.split()[8].split(';')

                        # use regular expressions to get the GeneID
                        geneID = re.search(r'GeneID:(\d+)', line)
                        if geneID:
                            value = geneID.group(1)
                            genomeRNAID[genomePath]['5S'] = value
      

        self.downloadrRNA(genomeRNAID)


        
    
    def getAccessions(self):
        ''' Get list of accessions from genomeInfo file '''

        fileName = f'gtrnadb_genomes.zip'
        #tmpAccess = f'genome_accessions_tmp.txt'
        tmpAccess = self.genomeInfo
        if os.path.exists(fileName):
            pass
        # else:
        
        #     accessions = self.genomeInfo['Rs_asm_acc'].tolist()
        #     with open(tmpAccess, 'w') as f:
        #         for acc in accessions:
        #             if acc != 'na':
        #                 f.write(acc.strip() + '\n')

        self.downloadNCBIgenomeInfo(tmpAccess, fileName, self.downloadDirectory)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python3 collect_genome_info.py genomeInfoFile downloadDirectory')
        sys.exit(1)
    #COLLECT_GENOME_INFO(sys.argv[1], sys.argv[2]).getAccessions()
    # get RNA
    COLLECT_GENOME_INFO(sys.argv[1], sys.argv[2]).getrRNA()
    #COLLECT_GENOME_INFO(sys.argv[1], sys.argv[2]).downloadNCBIgenomeInfo(accs='/projects/lowelab/users/jsleavit/git_repos/data/genome_accessions_tmp.txt',fileName='gtrnadb_genomes.zip', downloadDir=sys.argv[2])
