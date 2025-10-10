#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from math import log2


sns.set_theme(style="ticks")

plt.rcParams['axes.titlesize']     = 7        
plt.rcParams['axes.labelsize']     = 6        
plt.rcParams['xtick.labelsize']    = 6       
plt.rcParams['ytick.labelsize']    = 6       
plt.rcParams['legend.fontsize']    = 6       
plt.rcParams['figure.titlesize']   = 7   
plt.rcParams['font.size']          = 6
plt.rcParams['axes.linewidth']     = 0.6
plt.rcParams['xtick.major.width']  = 0.6
plt.rcParams['ytick.major.width']  = 0.6
plt.rcParams['xtick.minor.width']  = 0.4
plt.rcParams['ytick.minor.width']  = 0.4
plt.rcParams['xtick.major.size']   = 3
plt.rcParams['ytick.major.size']   = 3
plt.rcParams['xtick.minor.size']   = 2
plt.rcParams['ytick.minor.size']   = 2
plt.rcParams['legend.frameon']     = False
plt.rcParams["figure.figsize"]     = (3, 2)

red_shade = "#d62728"
grey_shade = "#7f7f7f"

custom_palette = sns.color_palette([red_shade, grey_shade], n_colors=2)
sns.set_palette(custom_palette)


def loadModules(tsv: str, min_size: int, max_size: int) -> list:
    modules = []
    with open(tsv, 'r') as modules_tsv:
        for line in modules_tsv:
            module = set(line.strip().split('\t'))
            if min_size <= len(module) <= max_size:
                modules.append(module)
    return modules


def sizeDistrib(modules: list, outFile: str, identifier: str):
    modulesSize = pd.DataFrame([len(m) for m in modules], columns=['Gene Set Size'])
    plt.clf()
    f, ax = plt.subplots()
    sns.despine(f)
    ax = sns.histplot(
        modulesSize,
        x="Gene Set Size",
        binwidth=5
    )
    ax.set_ylabel('Count (Gene Sets)')
    ax.set_title(f'{identifier} ({len(modules)} Gene Sets)')
    plt.tight_layout()
    plt.savefig(outFile)


def logSizeDistrib(modules: list, outFile: str, identifier: str):
    modulesSize = pd.DataFrame([log2(len(m)) for m in modules], columns=['log2(Gene Set Size)'])
    plt.clf()
    f, ax = plt.subplots()
    sns.despine(f)
    ax = sns.histplot(
        modulesSize,
        x="log2(Gene Set Size)",
        binwidth=0.5,
    )
    ax.set_ylabel('Count (Gene Sets)')
    ax.set(title=f'{identifier} ({len(modules)} Gene Sets)')
    plt.tight_layout()
    plt.savefig(outFile)


def pleiotropyDistrib(modules: list, outFile: str, identifier: str):
    genes = list(set.union(*modules))
    counts = [0 for i in range(len(genes))]
    i = -1
    for g in genes:
        i = i + 1
        for m in modules:
            if g in m:
                counts[i] = counts[i] + 1
    genesCount = pd.DataFrame(counts, index=genes, columns=["# Gene Occurrences in Sets"])
    plt.clf()
    f, ax = plt.subplots()
    sns.despine(f)
    ax = sns.histplot(
        genesCount,
        x="# Gene Occurrences in Sets",
        binwidth=10,
    )
    ax.set_ylabel('Count (Genes)')
    ax.set(title=f'{identifier} ({len(genes)} Genes)')
    ax.set_xlim(0, 350)
    plt.tight_layout()
    plt.savefig(outFile)


def logPleiotropyDistrib(modules: list, outFile: str, identifier: str):
    genes = list(set.union(*modules))
    counts = [0 for i in range(len(genes))]
    i = -1
    for g in genes:
        i = i + 1
        for m in modules:
            if g in m:
                counts[i] = counts[i] + 1
    counts = [log2(c) for c in counts if c > 0]
    genesCount = pd.DataFrame(counts, index=genes, columns=["log2(# Gene Occurrences in Sets)"])
    plt.clf()
    f, ax = plt.subplots()
    sns.despine(f)
    ax = sns.histplot(
        genesCount,
        x="log2(# Gene Occurrences in Sets)",
        binwidth=0.5,
    )
    ax.set_ylabel('Count (Genes)')
    ax.set(title=f'{identifier} ({len(genes)} Genes)')
    ax.set_xlim(0, 9)
    plt.tight_layout()
    plt.savefig(outFile)
    

def overlapRatioDistrib(modules: list, outFile: str, identifier: str):
    # for each module, compute the ratio between the amount of nodes shared with at least one other module, and its size
    overlapRatio = [
        round(len(
            modules[i].intersection(
                set.union(*[
                    x for x in modules if x != modules[i]
                ])
            )
        ) / len(modules[i]), 3) for i in range(len(modules))
    ]
    overlapRatio = pd.DataFrame(overlapRatio, columns=['Proportion of Shared Genes'])
    plt.clf()
    f, ax = plt.subplots()
    sns.despine(f)
    ax = sns.histplot(
        overlapRatio,
        x="Proportion of Shared Genes",
        binwidth=0.05
    )
    ax.set_ylabel('Count (Gene Sets)')
    ax.set(title=f'{identifier} ({len(modules)} Gene Sets)')
    plt.tight_layout()
    plt.savefig(outFile)

def jointSizeOverlapRatio(modules: list, outFile: str, identifier: str):
    df = pd.DataFrame(
        [len(m) for m in modules],
        columns=['Gene Set Size']
    )
    # for each module, compute the ratio between the amount of nodes shared with at least one other module, and its size
    overlapRatio = [
        round(len(
            modules[i].intersection(
                set.union(*[
                    x for x in modules if x != modules[i]
                ])
            )
        ) / len(modules[i]), 3) for i in range(len(modules))
    ]
    df['Proportion of Shared Genes'] = overlapRatio

    x = df['Gene Set Size']
    y = df['Proportion of Shared Genes']

    bw_x = 10
    bw_y = 0.05

    x_min = np.floor(x.min() / bw_x) * bw_x
    x_max = np.ceil(x.max() / bw_x) * bw_x
    y_min = np.floor(y.min() / bw_y) * bw_y
    y_max = np.ceil(y.max() / bw_y) * bw_y

    x_range = (x_min, x_max)
    y_range = (y_min, y_max)

    plt.clf()
    f, ax = plt.subplots()
    g = sns.jointplot(
        data=df,
        x='Gene Set Size',
        y='Proportion of Shared Genes',
        kind='hist',
        marginal_ticks=True,
        height=3,
        ratio=2,
        marginal_kws=dict(kde=False)
    )
    # Redraw marginals with their own bin widths but SAME ranges so edges align
    g.ax_marg_x.clear()
    sns.histplot(
        data=df, x="Gene Set Size",
        ax=g.ax_marg_x,
        binwidth=bw_x,
        binrange=x_range
    )
    g.ax_marg_y.clear()
    sns.histplot(
        data=df, y="Proportion of Shared Genes",
        ax=g.ax_marg_y,
        binwidth=bw_y,
        binrange=y_range
    )
    g.figure.suptitle(f"{identifier} ({len(modules)} Gene Sets)")
    g.ax_marg_y.set_ylabel("")
    g.ax_marg_x.set_xlabel("")
    plt.subplots_adjust(top=0.9)
    g.ax_marg_x.set_xlim(0, max(df['Gene Set Size']))
    g.ax_marg_y.set_ylim(0, 1)
    plt.tight_layout()
    plt.savefig(outFile)

def jointLogSizeOverlapRatio(modules: list, outFile: str, identifier: str):
    df = pd.DataFrame(
        [log2(len(m)) for m in modules],
        columns=['log(Gene Set Size)']
    )
    # for each modules, compute the ratio between the amount of nodes shared with at least one other module, and its size
    overlapRatio = [
        round(len(
            modules[i].intersection(
                set.union(*[
                    x for x in modules if x != modules[i]
                ])
            )
        ) / len(modules[i]), 3) for i in range(len(modules))
    ]
    df['Proportion of Shared Genes'] = overlapRatio

    x = df['log(Gene Set Size)']
    y = df['Proportion of Shared Genes']

    bw_x = 0.5
    bw_y = 0.05

    x_min = np.floor(x.min() / bw_x) * bw_x
    x_max = np.ceil(x.max() / bw_x) * bw_x
    y_min = np.floor(y.min() / bw_y) * bw_y
    y_max = np.ceil(y.max() / bw_y) * bw_y

    x_range = (x_min, x_max)
    y_range = (y_min, y_max)
    
    plt.clf()
    f, ax = plt.subplots()
    g = sns.jointplot(
        data=df,
        x='log(Gene Set Size)',
        y='Proportion of Shared Genes',
        kind='hist',
        marginal_ticks=True,
        height=3,
        ratio=2,
        marginal_kws=dict(kde=False),
        joint_kws=dict(
            binwidth=(bw_x, bw_y),
            binrange=[x_range, y_range],
        )
    )
    # Redraw marginals with their own bin widths but SAME ranges so edges align
    g.ax_marg_x.clear()
    sns.histplot(
        data=df, x="log(Gene Set Size)",
        ax=g.ax_marg_x,
        binwidth=bw_x,
        binrange=x_range
    )
    g.ax_marg_y.clear()
    sns.histplot(
        data=df, y="Proportion of Shared Genes",
        ax=g.ax_marg_y,
        binwidth=bw_y,
        binrange=y_range
    )
    g.figure.suptitle(f"{identifier} ({len(modules)} Gene Sets)")
    g.ax_marg_y.set_ylabel("")
    g.ax_marg_x.set_xlabel("")
    plt.subplots_adjust(top=0.9)
    g.ax_marg_x.set_xlim(0, max(df['log(Gene Set Size)']))
    g.ax_marg_y.set_ylim(0, 1)
    plt.tight_layout()
    plt.savefig(outFile)

def main():
    """
    modules_tsv_file = 'modules.tsv'
    sizeDistrib_pdf_file = 'sizeDistrib.pdf'
    logSizeDistrib_pdf_file = 'logSizeDistrib.pdf'
    pleiotropyDistrib_pdf_file = 'pleiotropyDistrib.pdf'
    overlapRatioDistrib_pdf_file = 'overlapRatioDistrib.pdf'
    jointSizeOverlapRatio_pdf_file = 'jointSizeOverlapRatio.pdf'
    jointLogSizeOverlapRatio_pdf_file = 'jointLogSizeOverlapRatio.pdf'
    """
    modules_tsv_file = sys.argv[1]
    logSizeDistrib_pdf_file = sys.argv[2]
    overlapRatioDistrib_pdf_file = sys.argv[3]
    pleiotropyDistrib_pdf_file = sys.argv[4]
    logPleiotropyDistrib_pdf_file = sys.argv[5]
    sizeDistrib_pdf_file = sys.argv[6]
    jointSizeOverlapRatio_pdf_file = sys.argv[7]
    jointLogSizeOverlapRatio_pdf_file = sys.argv[8]
    identifier = sys.argv[9]
    min_size = max(4, int(sys.argv[10]))
    max_size = min(512, int(sys.argv[11]))

    if identifier == 'wcsn_0.0_joined__wcmn_0.0_joined_merged':
        identifier = 'ProtVerse'
    elif identifier == 'wcsn_0.0_joined':
        identifier = 'Signalling Model'
    elif identifier == 'wcmn_0.0_joined':
        identifier = 'Metabolism Model'

    # read modules as list of sets
    modules = loadModules(tsv=modules_tsv_file, min_size=4, max_size=512)

    # plot histogram of module size
    #sizeDistrib(
    #    modules=modules,
    #    outFile=sizeDistrib_pdf_file,
    #    identifier=identifier
    #)
    #logSizeDistrib(
    #    modules=modules,
    #    outFile=logSizeDistrib_pdf_file,
    #    identifier=identifier
    #)

    # plot histogram of count of different modules per gene
    #pleiotropyDistrib(
    #    modules=modules,
    #    outFile=pleiotropyDistrib_pdf_file,
    #    identifier=identifier
    #)
    logPleiotropyDistrib(
        modules=modules,
        outFile=logPleiotropyDistrib_pdf_file,
        identifier=identifier
    )

    # plot histogram of proportion of shared nodes per module
    #overlapRatioDistrib(
    #    modules=modules,
    #    outFile=overlapRatioDistrib_pdf_file,
    #    identifier=identifier
    #)

    jointLogSizeOverlapRatio(
        modules=modules,
        outFile=jointLogSizeOverlapRatio_pdf_file,
        identifier=identifier
    )

    jointSizeOverlapRatio(
        modules=modules,
        outFile=jointSizeOverlapRatio_pdf_file,
        identifier=identifier
    )

if __name__ == '__main__':
    main()