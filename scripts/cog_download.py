#!/usr/bin/env python

import os
import json
import sys

class GET_COGs():
    def __init__(self, cog_json_path=''):
        if cog_json_path != '':
            self.cog_json = json.load(open(cog_json_path))

    def checkGtrnadbAssemblies(self, genome_accessions_list = []):
        ''' check if cog assembly is available in local ncbi database '''

        ''' 
        Note: 
        35 assemblies intersected in cog json
        35 assemblies intersected in gtRNAdb
        1309 assemblies in cog json
        4261 assemblies in genome path
        0.0078971119133574 percent of intersected genome assemblies found in genome path
        1.0127745595660463e-05 percent of intersected cog assemblies found in genome path
        '''

        cog_assemblies = [str(x['assembly']).strip('\n') for x in self.cog_json['cog-20.cog.csv']['features']]
        intersection = list(set(self.genome_accessions).intersection(set(cog_assemblies)))
        print(len(intersection), 'assemblies intersected in gtRNAdb')
        print(len(set(cog_assemblies)), 'assemblies in cog json')
        print(len(set(self.genome_accessions)), 'assemblies in genome path')
        print(len(intersection)/len(self.genome_accessions), 'percent of intersected genome assemblies found in genome path')
        print(len(intersection)/len(cog_assemblies), 'percent of intersected cog assemblies found in genome path')

    def download(self, cog, source = 'https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fasta/{0}.fa.gz', dlFile = '{0}.fasta.gz'):
        ''' download COG fasta files'''

        source = source.format(cog)
        dlFile = dlFile.format(cog)
        os.system('wget --output-document {0} {1}'.format(dlFile, source))
        return source, dlFile
    
    def downloadCOGs(self, outPath='', cogs = []):
        ''' download COGs in list '''

        downloaded_paths = []
        for cog in cogs:
            dl = f'{outPath}{cog}.fasta.gz'
            if not os.path.exists(dl):
                source, dl, = self.download(cog, dlFile = dl)
                os.system(f'gunzip --force {dl}')
            downloaded_paths.append(dl.replace('.gz', ''))
        return downloaded_paths
        

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('Usage: python cog_download.py download')
        print('Usage: python cog_download.py check')
        sys.exit(1)

    if sys.argv[1] == 'download':

        # download COGs
        print('Downloading COGs')
        out_path = '/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/fastas/'
        cogs_list = [x.strip('\n').split('\t')[0] for x in list(open('/projects/lowelab/users/jsleavit/git_repos/data/combined_trna_mod_database.tsv').readlines())][1:]
        downloadedCogs = GET_COGs().downloadCOGs(out_path, cogs_list)
    
    if sys.argv[1] == 'check':

        # check if cog assemblies are in local ncbi database
        cog_json_path='/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/cog-20.json'
        genome_accession_path='/projects/lowelab/users/jsleavit/git_repos/data/genome_accessions_tmp.txt'
        GET_COGs(cog_json_path).checkGtrnadbAssemblies(genome_accessions_list = [x.strip('\n') for x in list(open(genome_accession_path).readlines())])

