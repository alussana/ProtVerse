#!/usr/bin/env python3

# Example:
# cat gene_list.txt | gene2uniprot.py > gene2uniprot.txt

import urllib.parse
import urllib.request
import sys 

def gene2uniprot():
    nomenclature = 'ACC' # see https://www.uniprot.org/help/api_idmapping
    url = 'https://www.uniprot.org/uploadlists/'

    query = sys.stdin.read().replace('\n', ' ')

    params = {
        'from': nomenclature,
        'to': 'ENSEMBL_ID',
        #'to': 'CCDS_ID',
        'format': 'tab',
        'query': query,
        'taxon': '9606'
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
        print(response.decode('utf-8'))

if __name__ == '__main__':
    gene2uniprot()