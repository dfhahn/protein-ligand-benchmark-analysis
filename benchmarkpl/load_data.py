#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import re

import numpy as np
# import matplotlib
from matplotlib import pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from  plotly import colors
import pandas as pd
from rdkit import Chem

from rdkit.Chem import PandasTools
from IPython.core.display import HTML
from scipy.stats import norm

from PLBenchmarks import targets, ligands, edges

import pint
ureg = pint.UnitRegistry()

import benchmarkpl
path = benchmarkpl.__path__[0]
# In[ ]:


def getDetailedResults(target, forcefield='openff-1.0.0.offxml'):
    '''
    Get detailed results about pmx simulations with a specific target using forcefield.

    Parameters
    ----------
    target: string
    forcefield: string
    
    Returns
    -------
    pandas DataFrame with data
    '''
    # read in result file
    file_path = os.path.join(path, '..', '00_data', 'input', f'{target}_{forcefield}.csv')
    if not os.path.exists(file_path):
        return None
    res = pd.read_csv(file_path, header=[0,1,2], comment='#', skipinitialspace=True, sep=',')
    res.index = res.iloc[:,0].values
    
    res.rename(columns={'protein': 'complex'}, inplace=True)
    # read in exp. data
    edg = edges.EdgeSet(target)
    df = edg.get_dataframe(columns=[0,1, 'exp. DeltaG [kcal/mol]', 'exp. Error [kcal/mol]'])
    df.index = pd.Series(['edge_' + str(lig1) + '_' + str(lig2) for lig1, lig2 in zip(df[0].values, df[1].values)])
    # if target == 'jnk1':
    #     res.index = pd.Series([f'edge_{re.split("-|_", ind)[1]}_{re.split("-|_", ind)[3]}' for ind in res.index]) 
    #     df.index = pd.Series([f'edge_{re.split("-|_", ind)[1]}_{re.split("-|_", ind)[3]}' for ind in df.index])   
     
    res["lig1"] = df[0]
    res["lig2"] = df[1]
    # remove unit of calculated values
    res['exp_DDG'] = df['exp. DeltaG [kcal/mol]'].apply(lambda x: x.magnitude)
    res['dexp_DDG'] = df['exp. Error [kcal/mol]'].apply(lambda x: x.magnitude)

    return res


def getFFResults(target):
    vytas = pd.read_csv(f'/home/dhahn3/src/pmx/protLig_benchmark/ddg_data/{target}.dat', sep='\s+', header=None, comment='#',
                   names=['edge', 'exp', 'gaff', 'dgaff', 'cgenff', 'dcgenff', 'cons', 'dcons', 'fep5', 'dfep5', 'fep1', 'dfep1'])
    vytas['dexp']=pd.Series([0.0]*vytas.shape[0])
    newvytas = vytas.copy()
    
    newvytas.iloc[:,1:] = vytas.iloc[:,1:].apply(lambda x: pd.Series([round(ureg.Quantity(y, 'kJ / mole').to('kcal / mole').magnitude,2) for y in x]))
    newvytas.columns = pd.MultiIndex.from_arrays([np.array(vytas.columns), ['', 'exp'] + ['pmx'] * 6 + ['fep'] * 4 + ['exp'], [''] + ['kcal/mol'] * 12], names=['forcefield', 'method', 'unit'])
    
    # get own calculated values
    df = getResults(target)
    newvytas.index=['edge_'+ e for e in newvytas['edge']]
    newvytas.loc[:,('parsley', 'pmx', 'kcal/mol')] = df[('ddg_mean', '-', 'val')]
    newvytas.loc[:,('dparsley', 'pmx', 'kcal/mol')] = df[('ddg_mean', '-', 'err')]
    newvytas.loc[:,('exp', 'exp', 'kcal/mol')] = df['exp_DDG']
    newvytas.loc[:,('dexp', 'exp', 'kcal/mol')] = df['dexp_DDG']
    newvytas = newvytas.drop(columns=('edge', '', ''))
    newvytas.sort_index(axis=1, level=1, inplace=True, sort_remaining=False)
    return newvytas

