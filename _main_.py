# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 20:05:29 2020

@author: gosl241
"""

import sys,os
import argparse

import hyphalnet as hypnet

parser = argparse.ArgumentParser(description="""Get data from the proteomic /
                                 data commons and build community networks""")
                                 
def main():
    args = parser.parse_args()

if __name__ == '__main__':
    main()