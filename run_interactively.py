# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 11:29:06 2020
    Script to run code interactively at the terminal
@author: jgreensl
"""

import os # for local directory walk
import importlib # for importing at runtime
from utilities import fio

def selection(inlist):
    """
    list options, read selection, return selection
    """
    for i,item in enumerate(inlist):
        print("%d: %s"%(i,item))
    selection=str.lower(input("selection >> "))
    if str.isdigit(selection):
        selection=inlist[int(selection)]
    return selection

if __name__=='__main__':
    
    
    ## select from runnable scripts:
    scripts=[]
    for file in os.listdir("."):
        if file[-3:] == ".py":
            blacklist=['localtests','make_pft_files','run_interactively','tests']
            if file[:-3] not in blacklist:
                scripts.append(file[:-3])
    script=selection(scripts)
    
    ## import the selection
    print("INFO: importing ",script)
    importlib.import_module(script)
    help(script)
    ## print functions and docstring for selected script
    
    
    ## select from model runs
    model_runs=list(fio.model_outputs.keys())
    model_run=selection(model_runs)