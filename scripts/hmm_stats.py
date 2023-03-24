#!/usr/bin/env python

import sys
import os
import glob 

#hmmsearch_results = glob.glob('/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/hmms/*.hmmsearch.out')
hmmsearch_results = [x for x in glob.glob('/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/hmms/*.hmmsearch.out') if 'hmmsearch.out' in x]
for result in hmmsearch_results:
    cog = result.split('/')[-1].split('.')[0]
    if not os.path.exists('/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/hmms/'+cog):
        os.mkdir('/projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/hmmsearches/'+cog)
    os.system(f'mv {result} /projects/lowelab/users/jsleavit/git_repos/data/COG_ftp_files/hmmsearches/{cog}/')

