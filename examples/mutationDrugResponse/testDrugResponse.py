"""
EValuate the ability of hyphal networks to correlate with drug response
"""


import argparse
import sys
sys.path.insert(0, '../../')

from hyphalnet.hyha import hypha as hyp
import hyphalnet.mutations as mutes
import hyphalnet.hyphEnrich as hypex


parser = argparse.ArgumentParser(description="""Get data from mutations""")

def main():
    args = parser.parse_args()

    #get node-weighted data
    gfile='../../../OmicsIntegrator2/interactomes/inbiomap.9.12.2016.full.oi2'
    g = hyp.make_graph(gfile)
