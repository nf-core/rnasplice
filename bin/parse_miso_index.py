#!/usr/bin/env python
# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT

import glob
import shelve
import os
import argparse
import re

def main(prefix):

    # Check for non alphanumerics
    regex = re.compile('[@!#$%^&*()<>?/\|}{~:]')

    assert regex.search(prefix) == None, "Prefix contains special chars: " + prefix + ". Please remove special chars: [@!#$%^&*()<>?/\|}{~:]"

    # Open and extract shelve file contents
    genes_to_filenames = shelve.open(prefix + '/genes_to_filenames.shelve')
    data = [(gene, filename) for gene, filename in genes_to_filenames.iteritems()]
    genes_to_filenames.close()

    # Remove existing shelve files
    for filename in glob.glob(prefix + '/genes_to_filenames.shelve.*'):
        os.remove(filename)

    # Create new shelve file with relative paths
    genes_to_filenames = shelve.open(prefix + '/genes_to_filenames.shelve')

    for gene, filename in data:
        elements = filename.split(prefix)
        assert len(elements) == 2, "Filename contains" + (len(elements) - 1) + "instances of index delimiter. Expected only one."
        genes_to_filenames[gene] = './' + prefix + elements[1]

    genes_to_filenames.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Parse Miso index shelve files""")
    parser.add_argument(
        "-p",
        "--prefix",
        default="index",
        type=str,
        help="Folder where index shelve files stored",
    )

    args = parser.parse_args()
    main(args.prefix)
