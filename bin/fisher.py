#!/usr/bin/env python

# notes about scipy.stats.fisher_exact (FET)
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html

from scipy import stats
import sys
import pandas as pd
import numpy as np

def test_term(target: set, term: set, M: int, n_hypothesis: int):
    """
                             target set
                           - successes    - failures

    background - successes [  x           n - x     ]
               - failures  [N - x    M - (n + N) + x]
    """
    n = len(term)
    N = len(target)
    a = len(target.intersection(term))
    b = N - a
    c = n - a 
    #d = M - N - n + a
    d = M - a - b - c
    table = np.array([[a, b],
                      [c, d]])
    oddsr, p = stats.fisher_exact(table, alternative='greater')
    module_overlap_ratio = a / n
    target_overlap_ratio = a / N
    p_bonf = p * n_hypothesis
    return(round(oddsr, 5), round(p, 5), round(p_bonf, 5), a, b, c, d, round(module_overlap_ratio, 5), round(target_overlap_ratio, 5))

def main():
    """
    clusters_file = 'input/clusters.tsv'
    gene_set_file = 'input/genes.txt'
    gene_set_name = 'test'
    """
    clusters_file = sys.argv[1]
    gene_set_file = sys.argv[2]
    gene_set_name = sys.argv[3]
    
    # get gene set
    genes = []
    with open(gene_set_file) as gene_set:
        for line in gene_set:
            genes.append(line.strip())
    genes = set(genes)

    # get modules' names and modules' genes, as lists with consistent indeces
    terms_names = []
    terms_genes = []
    with open(clusters_file) as modules:
        for line in modules:
            try:
                term = line.strip().split('\t')
                terms_names.append(term[0])
                terms_genes.append(set(term[1:]))
            except:
                pass
    
    # union of the modules and the target genes
    u = set.union(*terms_genes)
    u = u.union(genes)
    M = len(u)
    
    # number of hypothesis to be tested for a given target gene set
    n_h = len(terms_genes)

    df_dict = {
        'set_name': [],
        'module_name': [],
        'oddsr': [],
        'module_fraction': [],
        'target_fraction': [],
        'fisher_p_greater': [],
        'fisher_p_greater_bonf': [],
        'a': [],
        'b': [],
        'c': [],
        'd': [],
        'genes_intersection': []
    }

    for i in range(n_h):
        term_name = terms_names[i]
        term_genes = terms_genes[i]
        oddsr, p, p_bonf, a, b, c, d, module_fraction, target_fraction = test_term(
            target=genes, term=term_genes, M=M, n_hypothesis=n_h)
        genes_intersect = term_genes.intersection(genes)
        df_dict['set_name'].append(gene_set_name)
        df_dict['module_name'].append(term_name)
        df_dict['oddsr'].append(oddsr)
        df_dict['module_fraction'].append(module_fraction)
        df_dict['target_fraction'].append(target_fraction)
        df_dict['fisher_p_greater'].append(p)
        df_dict['fisher_p_greater_bonf'].append(p_bonf)
        df_dict['a'].append(a)
        df_dict['b'].append(b)
        df_dict['c'].append(c)
        df_dict['d'].append(d)
        df_dict['genes_intersection'].append(genes_intersect)

    df = pd.DataFrame(df_dict).sort_values(by='fisher_p_greater_bonf')
    print(df.to_csv(sep='\t', header=True, index=False))

if __name__ == '__main__':
    """
    clusters_file = 'modules/wcsn_0.75/clusterone/clusters.tsv'
    gene_set_file = 'databases/koksal2018/diff_phos_genes.txt'
    """
    main()
