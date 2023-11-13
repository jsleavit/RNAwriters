#!/bin/usr/env python

import sys
import os
import pandas as pd
import gzip
import glob
import io


class FastAreader():
    def __init__(self, fname=''):
        ''' contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        header = ''
        sequence = ''

        with self.doOpen() as fileH:
            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


# COG BLAST
################################################################################
class COG_BLAST():
    def __init__(self,  genomeDir =  '/projects/lowelab/users/jsleavit/git_repos/data/ncbi/ncbi_dataset/data/',
                  outPath = '/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/',
                  blastReference = '/projects/lowelab/users/jsleavit/git_repos/data/uniprot_trna_mod_enzymes_references_query_set.fasta',
                    overwrite=False):

        self.genomeInfo = pd.read_csv('/projects/lowelab/users/jsleavit/git_repos/data/Gtrnadb_prokaryotic_info_with_full_taxonomy_updated_asm_accessions.tsv', sep = '\t', index_col = 0)
        self.out_path_root = outPath
        self.overwrite = overwrite
        self.blast_reference = blastReference
        self.genome_dir = genomeDir
    
    def makeDirectory(self, directory):
        ''' make a directory if it does not exist '''
        
        if not os.path.exists(directory):
            os.makedirs(directory)
        if os.path.exists(directory):
            return True
        
    def makeBlastDB(self):
        ''' make blast db from reference fasta '''

        for dir in ['blast_db/', 'arch_hmm_blast_results/']:
            if self.makeDirectory(self.out_path_root + dir):
                pass
            else:
                self.makeDirectory(self.out_path_root + dir)

        if os.path.exists(self.out_path_root + 'blast_db/uniprot_trna_mod_enzymes_references_query_set.fasta.nhr') and not self.overwrite:
            print('blast db exists')
            pass
        else:
            os.system(f'makeblastdb -in {self.blast_reference} -dbtype prot -out {self.out_path_root}blast_db/uniprot_trna_mod_enzymes_references_query_set.fasta')
    
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
    
    def mapUniprot(self):
        ''' build a dictionary of uniprot ids to gene names '''

        self.uniprot_map = {}
        seqs = self.readFasta(self.blast_reference)
        for seq in seqs.keys():
            self.uniprot_map[seq.split('|')[1]] = (seq.split('|')[2].split()[0], ' '.join([x for x in seq.split()[2:] if '=' not in x]))
        return self.uniprot_map
        
    
    def removeTempFiles(self):
        ''' remove temporary files '''

        os.system(f'rm -r {self.out_path_root}arch_temp/')
                    
    def blastHMMhits(self, highConfHitCutoff = 1e-40, eValueCutoff = 1e-3, numThreads = 20):
        ''' blast hmmsearch hits against reference fasta '''

        blast_db = self.out_path_root + 'blast_db/uniprot_trna_mod_enzymes_references_query_set.fasta'

        hmmsearch_stat_dir_list = glob.glob(self.out_path_root + 'hmmsearch_stats/*')
    
        hits = pd.DataFrame()

        uniprot_ids = self.mapUniprot()
        subset = {'Haloferax volcanii DS2': 'GCA_000025685.1', 'Methanocaldococcus jannaschii DSM 2661': 'GCA_000091665.1', 'Methanococcus maripaludis C5': 'GCA_000016125.1', 'Methanococcus maripaludis C6': 'GCA_000018485.1', 'Methanococcus maripaludis C7': 'GCA_000017225.1', 'Methanococcus maripaludis S2': 'GCA_000011585.1', 'Methanococcus maripaludis X1': 'GCA_000220645.1', 'Pyrococcus furiosus COM1': 'GCA_000275605.1', 'Pyrococcus furiosus DSM 3638': 'GCA_000007305.1', 'Sulfolobus acidocaldarius DSM 639': 'GCA_000012285.1', 'Sulfolobus acidocaldarius N8': 'GCA_000340315.1', 'Sulfolobus acidocaldarius Ron12/I': 'GCA_000338775.1', 'Sulfolobus acidocaldarius SUSAZ': 'GCA_000508305.1', 'Thermococcus kodakarensis KOD1': 'GCA_000009965.1'}
        subset_accs = [x.strip() for x in subset.values()]
        
        # open each directory in hmmsearch_stats and each file in that directory and read in the gene name to blast
        #for hmmsearch_dir in os.listdir(self.out_path_root + 'hmmsearch_stats/'):
        for hmmsearch_dir in hmmsearch_stat_dir_list:

            acc = hmmsearch_dir.split('/')[-1].strip()
            

            if acc in subset_accs:

                hmmsearch_dir = hmmsearch_dir.split('/')[-1].strip()
            

                # # check if the hmmsearch_dir is a directory
                # if os.path.isdir(self.out_path_root + 'hmmsearch_stats/' + hmmsearch_dir):

                # open each file in hmmsearch_dir
                for hmmsearch_file in os.listdir(self.out_path_root + 'hmmsearch_stats/' + hmmsearch_dir):

                    with open(self.out_path_root + 'hmmsearch_stats/' + hmmsearch_dir + '/' + hmmsearch_file, 'r') as f:
                        for line in f:
                            cog = hmmsearch_file.split('.')[3]

                            # clean up the line
                            gene, e_val, product, phylum, family, species = [str(x).strip() for x in line.split('\t')[3:]]

                            # make a temporary directory for blasting
                            if self.makeDirectory(self.out_path_root + 'arch_temp/'):
                                pass
                            else:
                                self.makeDirectory(self.out_path_root + 'arch_temp/')
                            
                            # get all the protein sequences from the each cog in the hmmsearch_dir
                            with open(self.out_path_root + f'arch_temp/{hmmsearch_dir}.fasta', 'a') as big_fasta_out:
                                FA = FastAreader(self.genome_dir + hmmsearch_dir + '/protein.faa')
                                for head, seq in FA.readFasta():
                                    if gene in head.split()[0]:
                                        #gene_id = head.split()[0]
                                        annotated_product = '_'.join(head.split()[1:])
                                        #gene_with_out_spaces = '_'.join(head.split())
                                        big_fasta_out.write('>{}|{}|{}|{}|{}|{}\n{}\n'.format(gene,cog,e_val,hmmsearch_dir,phylum,annotated_product,seq))
                            query_path = self.out_path_root + f'arch_temp/{hmmsearch_dir}.fasta'

                # blast gene against reference fasta
                print(f'blasting {hmmsearch_dir} against {blast_db}')
                os.system(f'blastp -query {query_path} -db {blast_db} -num_threads {numThreads} -outfmt 6 -out {self.out_path_root}arch_hmm_blast_results/{hmmsearch_dir}.blast.out')
                
                # remove the temporary directory
                self.removeTempFiles()

                # build the blast results database
                gene_hits = pd.read_csv(self.out_path_root + f'arch_hmm_blast_results/{hmmsearch_dir}.blast.out', comment = '#', sep = '\t', header = None).rename(columns = {0: 'query', 1: 'target', 2: 'percent_identity',  10: 'e_value', 11: 'bit_score'})
                gene_hits = gene_hits[gene_hits['e_value'] <= eValueCutoff]
                gene_hits['target_product'] = gene_hits['target'].apply(lambda x: uniprot_ids[x][1])
                hits = pd.concat([hits, gene_hits], ignore_index = True)

        # write out the blast results
        hits.to_csv(self.out_path_root + 'arch_hmm_blast_results/hmm_blast_results.tsv', index = False, header = True, sep = '\t')
    
if __name__ == '__main__':
    COG_BLAST().blastHMMhits()
   