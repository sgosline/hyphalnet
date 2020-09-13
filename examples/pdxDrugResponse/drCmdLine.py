"""
Evaluate the ability of hyphal networks to correlate with drug response
Uses command line version of hyphalnet
"""


import argparse
import os
#sys.path.insert(0, '../../')


import pandas as pd
import synapseclient as synapse

parser = argparse.ArgumentParser(description="""Get data from mutations""")


def getMutationalData():
    syn = synapse.login()
    tab = syn.tableQuery("SELECT * FROM syn22266278").asDataFrame()
    tab = tab.rename(columns={'specimenID':'Sample','AD':'Value'})
    return tab

def main():
    args = parser.parse_args()
    beta = 0.5
    #get mutational data
    mvals = getMutationalData().to_csv("mutData.csv")

    key = 'mpnstPDXmuts'
    os.system('python3 ../../bin/createHyphaFromProts.py --output '+key+' --inputData mutData.csv --graph ../../data/igraphPPI.pkl')

if __name__=='__main__':
    main()
