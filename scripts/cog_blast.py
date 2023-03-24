#!/bin/usr/env python

import sys
import os
import pandas as pd
import gzip
import glob
import io

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

        for dir in ['blast_db/', 'hmm_blast_results/']:
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
            f.close

    def getProteinSequence(self, genome, gene):
        ''' get protein sequence from genome fasta file and write a temporary fasta file for blasting '''

        # make a temporary directory for blasting
        if self.makeDirectory(self.out_path_root + 'temp/'):
            pass
        else:
            self.makeDirectory(self.out_path_root + 'temp/')

        genome_fasta = self.genome_dir + genome + '/protein.faa'
        seqs = self.readFasta(genome_fasta)
        
        # check only the first element of each key in the dictionary to see if it matches the gene name
        if gene in map(lambda x: x.split()[0], seqs.keys()):
            gene = list(filter(lambda x: x.split()[0] == gene, seqs.keys()))[0]
            gene_id = gene.split()[0]
            product = ' '.join(gene.split()[1:])
            gene_with_out_spaces = '_'.join(gene.split())
            self.writeFasta({gene_with_out_spaces: seqs[gene]}, self.out_path_root + f'temp/{genome}.{gene_id}.fasta')
        return self.out_path_root + f'temp/{genome}.{gene_id}.fasta', product
    
    def mapUniprot(self):
        ''' build a dictionary of uniprot ids to gene names '''

        self.uniprot_map = {}
        seqs = self.readFasta(self.blast_reference)
        for seq in seqs.keys():
            self.uniprot_map[seq.split('|')[1]] = (seq.split('|')[2].split()[0], ' '.join([x for x in seq.split()[2:] if '=' not in x]))
        return self.uniprot_map
        
    
    def removeTempFiles(self):
        ''' remove temporary files '''

        os.system(f'rm -r {self.out_path_root}temp/')
                    
    def blastHMMhits(self, highConfHitCutoff = 1e-40, eValueCutoff = 1e-3, numThreads = 8):
        ''' blast hmmsearch hits against reference fasta '''

        blast_db = self.out_path_root + 'blast_db/uniprot_trna_mod_enzymes_references_query_set.fasta'

        hits = pd.DataFrame()

        uniprot_ids = self.mapUniprot()

        # open each directory in hmmsearch_stats and each file in that directory and read in the gene name to blast
        for hmmsearch_dir in os.listdir(self.out_path_root + 'hmmsearch_stats/'):

            # check if the hmmsearch_dir is a directory
            if os.path.isdir(self.out_path_root + 'hmmsearch_stats/' + hmmsearch_dir):

                # make a new directory for each hmmsearch_dir in hmm_blast_results
                if self.makeDirectory(self.out_path_root + f'hmm_blast_results/{hmmsearch_dir}/'):
                    pass
                else:
                    self.makeDirectory(self.out_path_root + f'hmm_blast_results/{hmmsearch_dir}/')

                # open each file in hmmsearch_dir
                for hmmsearch_file in os.listdir(self.out_path_root + 'hmmsearch_stats/' + hmmsearch_dir):

                    with open(self.out_path_root + 'hmmsearch_stats/' + hmmsearch_dir + '/' + hmmsearch_file, 'r') as f:
                        for line in f:
                            cog = hmmsearch_file.split('.')[3]

                            # make a new directory for each cog in hmm_blast_results
                            if self.makeDirectory(self.out_path_root + f'hmm_blast_results/{hmmsearch_dir}/{cog}/'):
                                pass
                            else:
                                self.makeDirectory(self.out_path_root + f'hmm_blast_results/{hmmsearch_dir}/{cog}/')

                            # clean up the line
                            gene, e_val, product, phylum, family, species = [str(x).strip() for x in line.split('\t')[3:]]

                            # get the protein sequence from the genome fasta file
                            gene_path, product = self.getProteinSequence(hmmsearch_dir, gene)

                            # blast gene against reference fasta
                            os.system(f'blastp -query {gene_path} -db {blast_db} -num_threads {numThreads} -outfmt 6 -out {self.out_path_root}hmm_blast_results/{hmmsearch_dir}/{cog}/{gene}.blast.out')

                            # only read to csv if the file is not empty
                            if open(self.out_path_root + f'hmm_blast_results/{hmmsearch_dir}/{cog}/{gene}.blast.out', 'r').read() != '':
                                gene_hits = pd.read_csv(self.out_path_root + f'hmm_blast_results/{hmmsearch_dir}/{cog}/{gene}.blast.out', comment = '#', sep = '\t', header = None).rename(columns = {0: 'query', 1: 'target', 2: 'percent_identity',  10: 'e_value', 11: 'bit_score'})
                                
                                gene_hits = gene_hits[gene_hits['e_value'] <= eValueCutoff]

                                gene_hits['query'] = gene
                                gene_hits['target_source'] = gene_hits['target'].apply(lambda x: uniprot_ids[x][0])
                                gene_hits['cog'] = cog
                                gene_hits['phylum'] = phylum
                                gene_hits['family'] = family
                                gene_hits['species'] = species
                                gene_hits['product'] = product
                                gene_hits['target_product'] = gene_hits['target'].apply(lambda x: uniprot_ids[x][1])
                                gene_hits['hmm_eval'] = e_val

                                hits = pd.concat([hits, gene_hits])


                # clear the temporary directory for next species 
                self.removeTempFiles()

            # if it is not a directory, continue
            else:
                continue

        # write the dataframe out to a tsv file]
        try: 
            print('writing to pickle')
            hits.to_pickle(self.out_path_root + 'hmm_blast_results/hmm_blast_results.pkl')
        except:
            print('writing to tsv')
            hits.to_csv(self.out_path_root + 'hmm_blast_results/hmm_blast_results.tsv', sep = '\t', index = False)


if __name__ == '__main__':
    COG_BLAST().blastHMMhits()
   