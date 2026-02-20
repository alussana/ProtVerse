#!/usr/bin/env python3

import sys
import numpy as np


def main():
	in_path = sys.argv[1]

	data = np.atleast_2d(np.genfromtxt(in_path, delimiter="\t", dtype=str))
	n1, n2 = data[:, 0], data[:, 1]
	rest = data[:, 2:]

	mask = n1 <= n2
	out = np.column_stack([np.where(mask, n1, n2),
                           np.where(mask, n2, n1),
                           rest])

	np.savetxt(sys.stdout, out, delimiter="\t", fmt="%s")


if __name__ == '__main__':
	main()
