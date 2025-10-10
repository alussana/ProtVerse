#!/usr/bin/env python3

# Example:
# cat gene_list.txt | gene2uniprot.py <FROM_NOMENCLATURE> <TO_NOMENCLATURE> <TAXONOMY_ID> > dict.txt

# see https://www.uniprot.org/help/api_idmapping for nomenclatures and taxonomy ids

import urllib.parse
import urllib.request
import sys 

def gene2uniprot():
    from_ = sys.argv[1]
    to_= sys.argv[2]
    taxon = sys.argv[3]
    url = 'https://www.uniprot.org/uploadlists/'

    query = sys.stdin.read().replace('\n', ' ')

    params = {
        'from': from_,
        'to': to_,
        'format': 'tab',
        'query': query,
        'taxon': taxon
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
        print(response.decode('utf-8'))

if __name__ == '__main__':
    gene2uniprot()