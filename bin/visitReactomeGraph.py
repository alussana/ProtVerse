#!/usr/bin/env python

import sys

def visitNode(edges: list, node: str):
    child_nodes = []
    
    for i in range(len(edges)):
        if edges[i][0] == node:
            child_nodes.append(edges[i][1])
    
    if len(child_nodes) == 0:
        print(node)
    else:
        for i in range(len(child_nodes)):
            visitNode(edges, child_nodes[i])

def main():
    file = sys.argv[1]
    node = sys.argv[2]
    
    graph = []
    
    with open(file) as fh:
        for line in fh:
            edge = line.strip().split('\t')
            graph.append(edge)
    
    visitNode(graph, node)

if __name__ == '__main__':
    main()