import os
import numpy as np
from scipy.stats import norm
import pandas as pd

import plotly.graph_objects as go
from plotly.subplots import make_subplots
from  plotly import colors
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor, rdFMCS, TemplateAlign
from rdkit.Chem.Draw import rdMolDraw2D
from svgutils import transform as sg
from PLBenchmarks import targets, ligands, edges


newcolors = []
c = [colors.hex_to_rgb(c) for c in colors.cyclical.mrybm]
#c = [colors.unlabel_rgb(c) for c in colors.sequential.Bluyl]

for i, e in enumerate(c[8:14]):
    newcolors.append(e)
    newcolors.append(colors.find_intermediate_color(e, c[i+9], 0.5))
newcolors = [ colors.unconvert_from_RGB_255(c) for c in newcolors ]
newcolorsH = []
# c = [colors.unlabel_rgb(c) for c in colors.sequential.Peach]

for i, e in enumerate(c[:6]):
    newcolorsH.append(e)
    newcolorsH.append(colors.find_intermediate_color(e, c[i+1], 0.5))
newcolorsH = [ colors.unconvert_from_RGB_255(c) for c in newcolorsH ]

def drawPerturbation(m1, m2, pairs, target='', n1='', n2='', text=''):
    d2d = rdMolDraw2D.MolDraw2DSVG(600,300,300,300)
    drawoptions = d2d.drawOptions()
    drawoptions.padding = 0.15
    drawoptions.prepareMolsBeforeDrawing = True
#     drawoptions.addAtomIndices = True
#     drawoptions.atomsLabels = True
#     d2d.SetDrawOptions(drawoptions)
#     drawoptions.centreMoleculesBeforeDrawing = False
    # Generate 2D Depictions
    mcs = rdFMCS.FindMCS([m1,m2],  maximizeBonds=True, bondCompare=rdFMCS.BondCompare.CompareAny, atomCompare=rdFMCS.AtomCompare.CompareElements,
                         timeout=10) # find maximum common substructure (MCS)
    pattern = Chem.MolFromSmarts(mcs.smartsString) # generate mol with MCS
    match = m1.GetSubstructMatch(pattern) # get substruct match of MCS with first mol
    conf = m1.GetConformer(0) # get first molecule's conformer

    cnf = Chem.Conformer(len(match)) # generate a conformer for pattern
    [ cnf.SetAtomPosition(i, conf.GetAtomPosition(j)) for i, j in enumerate(match) ] # set pattern's atom position to first mol's atom positions
    pattern.AddConformer(cnf) # add conformer to pattern
    
    Chem.rdDepictor.GenerateDepictionMatching3DStructure(pattern, pattern) # generate depict of pattern matching 3D struct
    TemplateAlign.AlignMolToTemplate2D(m1,pattern, clearConfs=True) # generate depict for m1 aligned to pattern
    TemplateAlign.AlignMolToTemplate2D(m2,pattern, clearConfs=True) # generate depict for m2 aligned to pattern

    specialHighlight = []
    for i, p in enumerate(pairs):
        a1 = m1.GetAtomWithIdx(int(p[0]))
        a2 = m2.GetAtomWithIdx(int(p[1]))
        sameElement = a1.GetSymbol() == a2.GetSymbol()
        bothAromatic = a1.GetIsAromatic() and a2.GetIsAromatic()
        bothNotAromatic = not a1.GetIsAromatic() and not a2.GetIsAromatic()
        sameHybridization = a1.GetHybridization() == a2.GetHybridization()
        specialHighlight.append(not sameElement or not sameHybridization or not (bothAromatic or bothNotAromatic))
    
    highlightAtoms = [[int(p[0]) for p in pairs],
                      [int(p[1]) for p in pairs]]
    highlightAtomColors = [{int(p[0]): newcolorsH[i%len(newcolorsH)] if specialHighlight[i]  else newcolors[i%len(newcolors)] for i, p in enumerate(pairs)},
                           {int(p[1]): newcolorsH[i%len(newcolorsH)] if specialHighlight[i]  else newcolors[i%len(newcolors)] for i, p in enumerate(pairs)}]
    highlightAtomRadii = [{int(p[0]): .5 if specialHighlight[i]  else .33 for i, p in enumerate(pairs)},
                           {int(p[1]): .5 if specialHighlight[i]  else .33 for i, p in enumerate(pairs)}]
    
#     for atom in m1.GetAtoms():
#         atom.SetProp('atomLabel',str(atom.GetIdx()))
#     for atom in m2.GetAtoms():
#         atom.SetProp('atomLabel',str(atom.GetIdx()))        
    d2d.DrawMolecules([m1, m2],  
                      highlightAtoms=highlightAtoms, 
                      highlightAtomColors=highlightAtomColors, 
                      highlightAtomRadii=highlightAtomRadii, 
                      highlightBonds=[[],[]])

#     conf = m2.GetConformer(0)
#     crds=np.array([list(conf.GetAtomPosition(i)) for i in range(m2.GetNumAtoms())])
#     xrange = np.max(crds[:,0]) - np.min(crds[:,0])
#     yrange = np.max(crds[:,1]) - np.min(crds[:,1])
#     textx = -xrange*0.7
#     texty = yrange*0.7
# #d2d.SetOffset(0,0)
#     drawoptions.setSymbolColour((1,1,1))
#     d2d.SetFontSize(d2d.FontSize()*1.5)
#     d2d.DrawString(f'{target}', Point2D(textx, texty))
#     d2d.DrawString(f'{n1}->{n2}', Point2D(textx, texty-.75))
#     print(list(d2d.Offset()))
#     print(list(d2d.GetDrawCoords(Point2D(np.min(crds[:,0]), np.max(crds[:,1])))))
    d2d.FinishDrawing()
    s = d2d.GetDrawingText()
#     s = s.replace('svg:','')
    fig = sg.fromstring(s)
    
    label = sg.TextElement(300, 20, f'{target} {n1} -> {n2} ', size=20, 
                       font='sans-serif', anchor='middle', color='black')
    label1 = sg.TextElement(150, 285, f'{n1}', size=15, 
                       font='sans-serif', anchor='middle', color='black')
    label2 = sg.TextElement(450, 285, f'{n2}', size=15, 
                       font='sans-serif', anchor='middle', color='black')
    fig.append(label) 
    fig.append(label1) 
    fig.append(label2)
    return fig.to_str().decode("utf-8") 

def drawPerturbationInverted(m1, m2, pairs, target='', n1='', n2='', text=''):
    d2d = rdMolDraw2D.MolDraw2DSVG(600,300,300,300)
    drawoptions = d2d.drawOptions()
    drawoptions.padding = 0.15
    drawoptions.prepareMolsBeforeDrawing = True
#     drawoptions.addAtomIndices = True
#     drawoptions.atomsLabels = True
#     d2d.SetDrawOptions(drawoptions)
#     drawoptions.centreMoleculesBeforeDrawing = False
    # Generate 2D Depictions
    mcs = rdFMCS.FindMCS([m1,m2],  maximizeBonds=True, bondCompare=rdFMCS.BondCompare.CompareAny, atomCompare=rdFMCS.AtomCompare.CompareElements,
                         timeout=10) # find maximum common substructure (MCS)
    pattern = Chem.MolFromSmarts(mcs.smartsString) # generate mol with MCS
    match = m1.GetSubstructMatch(pattern) # get substruct match of MCS with first mol
    conf = m1.GetConformer(0) # get first molecule's conformer

    cnf = Chem.Conformer(len(match)) # generate a conformer for pattern
    [ cnf.SetAtomPosition(i, conf.GetAtomPosition(j)) for i, j in enumerate(match) ] # set pattern's atom position to first mol's atom positions
    pattern.AddConformer(cnf) # add conformer to pattern
    
    Chem.rdDepictor.GenerateDepictionMatching3DStructure(pattern, pattern) # generate depict of pattern matching 3D struct
    TemplateAlign.AlignMolToTemplate2D(m1,pattern, clearConfs=True) # generate depict for m1 aligned to pattern
    TemplateAlign.AlignMolToTemplate2D(m2,pattern, clearConfs=True) # generate depict for m2 aligned to pattern

    specialHighlight = []
    for i, p in enumerate(pairs):
        a1 = m1.GetAtomWithIdx(int(p[0]))
        a2 = m2.GetAtomWithIdx(int(p[1]))
        sameElement = a1.GetSymbol() == a2.GetSymbol()
        bothAromatic = a1.GetIsAromatic() and a2.GetIsAromatic()
        bothNotAromatic = not a1.GetIsAromatic() and not a2.GetIsAromatic()
        sameHybridization = a1.GetHybridization() == a2.GetHybridization()
        specialHighlight.append(not sameElement or not sameHybridization or not (bothAromatic or bothNotAromatic))
    highlightAtoms = [[int(p[0])  for i, p in enumerate(pairs) if specialHighlight[i]],
                      [int(p[1])  for i, p in enumerate(pairs) if specialHighlight[i]]]
    highlightAtomColors = [{int(p[0]): newcolors[i%len(newcolors)]for i, p in enumerate(pairs) if specialHighlight[i]},
                           {int(p[1]): newcolors[i%len(newcolors)] for i, p in enumerate(pairs) if specialHighlight[i]}]
    highlightAtomRadii = [{int(p[0]): .33 for i, p in enumerate(pairs) if specialHighlight[i]},
                           {int(p[1]): .33 for i, p in enumerate(pairs) if specialHighlight[i]}]

    num1 = m1.GetNumAtoms()
    for i in range(num1):
        if i not in pairs[:,0]:
            highlightAtoms[0].append(i)
            highlightAtomColors[0][i]=newcolorsH[5]
            highlightAtomRadii[0][i]=.33
    num2 = m2.GetNumAtoms()
    for i in range(num2):
        if i not in pairs[:,1]:
            highlightAtoms[1].append(i)
            highlightAtomColors[1][i]=newcolorsH[5]
            highlightAtomRadii[1][i]=.33
    
#     for atom in m1.GetAtoms():
#         atom.SetProp('atomLabel',str(atom.GetIdx()))
#     for atom in m2.GetAtoms():
#         atom.SetProp('atomLabel',str(atom.GetIdx()))        
    d2d.DrawMolecules([m1, m2],  
                      highlightAtoms=highlightAtoms, 
                      highlightAtomColors=highlightAtomColors, 
                      highlightAtomRadii=highlightAtomRadii, 
                      highlightBonds=[[],[]])

#     conf = m2.GetConformer(0)
#     crds=np.array([list(conf.GetAtomPosition(i)) for i in range(m2.GetNumAtoms())])
#     xrange = np.max(crds[:,0]) - np.min(crds[:,0])
#     yrange = np.max(crds[:,1]) - np.min(crds[:,1])
#     textx = -xrange*0.7
#     texty = yrange*0.7
# #d2d.SetOffset(0,0)
#     drawoptions.setSymbolColour((1,1,1))
#     d2d.SetFontSize(d2d.FontSize()*1.5)
#     d2d.DrawString(f'{target}', Point2D(textx, texty))
#     d2d.DrawString(f'{n1}->{n2}', Point2D(textx, texty-.75))
#     print(list(d2d.Offset()))
#     print(list(d2d.GetDrawCoords(Point2D(np.min(crds[:,0]), np.max(crds[:,1])))))
    d2d.FinishDrawing()
    s = d2d.GetDrawingText()
#     s = s.replace('svg:','')
    fig = sg.fromstring(s)
    
    label = sg.TextElement(300, 20, f'{target} {n1} -> {n2} ', size=20, 
                       font='sans-serif', anchor='middle', color='black')
    label1 = sg.TextElement(150, 285, f'{n1}', size=15, 
                       font='sans-serif', anchor='middle', color='black')
    label2 = sg.TextElement(450, 285, f'{n2}', size=15, 
                       font='sans-serif', anchor='middle', color='black')
    fig.append(label) 
    fig.append(label1) 
    fig.append(label2)
    return fig.to_str().decode("utf-8") 

def add_test_to_perturbation(fig_string, text):
    fig = sg.fromstring(s)
    
    label = sg.TextElement(300, 35, f'{text}', size=15, 
                       font='sans-serif', anchor='middle', color='black')
    fig.append(label) 
    return fig.to_str().decode("utf-8") 

def drawPerturbationBare(m1, m2, pairs, target='', n1='', n2='', text=''):
    d2d = rdMolDraw2D.MolDraw2DSVG(600,300,300,300)
    drawoptions = d2d.drawOptions()
    drawoptions.padding = 0.15
    drawoptions.prepareMolsBeforeDrawing = True
#     drawoptions.addAtomIndices = True
#     drawoptions.atomsLabels = True
#     d2d.SetDrawOptions(drawoptions)
#     drawoptions.centreMoleculesBeforeDrawing = False
    # Generate 2D Depictions
#         atom.SetProp('atomLabel',str(atom.GetIdx()))
#     for atom in m2.GetAtoms():
#         atom.SetProp('atomLabel',str(atom.GetIdx()))        
    Chem.rdDepictor.Compute2DCoords(m1)
    Chem.rdDepictor.Compute2DCoords(m2)
    d2d.DrawMolecules([m1, m2])

#     conf = m2.GetConformer(0)
#     crds=np.array([list(conf.GetAtomPosition(i)) for i in range(m2.GetNumAtoms())])
#     xrange = np.max(crds[:,0]) - np.min(crds[:,0])
#     yrange = np.max(crds[:,1]) - np.min(crds[:,1])
#     textx = -xrange*0.7
#     texty = yrange*0.7
# #d2d.SetOffset(0,0)
#     drawoptions.setSymbolColour((1,1,1))
#     d2d.SetFontSize(d2d.FontSize()*1.5)
#     d2d.DrawString(f'{target}', Point2D(textx, texty))
#     d2d.DrawString(f'{n1}->{n2}', Point2D(textx, texty-.75))
#     print(list(d2d.Offset()))
#     print(list(d2d.GetDrawCoords(Point2D(np.min(crds[:,0]), np.max(crds[:,1])))))
    d2d.FinishDrawing()
    s = d2d.GetDrawingText()
#     s = s.replace('svg:','')
    fig = sg.fromstring(s)
    
    label = sg.TextElement(300, 20, f'{target}: {n1} -> {n2}', size=20, 
                       font='sans-serif', anchor='middle', color='black')
    label2 = sg.TextElement(300, 35, f'{text}', size=15, 
                       font='sans-serif', anchor='middle', color='black')
    label3 = sg.TextElement(150, 285, f'{n1}', size=15, 
                       font='sans-serif', anchor='middle', color='black')
    label4 = sg.TextElement(450, 285, f'{n2}', size=15, 
                       font='sans-serif', anchor='middle', color='black')
    fig.append(label) 
    fig.append(label2) 
    fig.append(label3)
    fig.append(label4)
    return fig.to_str().decode("utf-8") 

def drawPerturbationBareInverted(m1, m2, pairs, target='', n1='', n2='', text=''):
    d2d = rdMolDraw2D.MolDraw2DSVG(600,300,300,300)
    drawoptions = d2d.drawOptions()
    drawoptions.padding = 0.15
    drawoptions.prepareMolsBeforeDrawing = True
#     drawoptions.addAtomIndices = True
#     drawoptions.atomsLabels = True
#     d2d.SetDrawOptions(drawoptions)
#     drawoptions.centreMoleculesBeforeDrawing = False
    # Generate 2D Depictions
#         atom.SetProp('atomLabel',str(atom.GetIdx()))
#     for atom in m2.GetAtoms():
#         atom.SetProp('atomLabel',str(atom.GetIdx()))        
    Chem.rdDepictor.Compute2DCoords(m1)
    Chem.rdDepictor.Compute2DCoords(m2)
    d2d.DrawMolecules([m1, m2])

#     conf = m2.GetConformer(0)
#     crds=np.array([list(conf.GetAtomPosition(i)) for i in range(m2.GetNumAtoms())])
#     xrange = np.max(crds[:,0]) - np.min(crds[:,0])
#     yrange = np.max(crds[:,1]) - np.min(crds[:,1])
#     textx = -xrange*0.7
#     texty = yrange*0.7
# #d2d.SetOffset(0,0)
#     drawoptions.setSymbolColour((1,1,1))
#     d2d.SetFontSize(d2d.FontSize()*1.5)
#     d2d.DrawString(f'{target}', Point2D(textx, texty))
#     d2d.DrawString(f'{n1}->{n2}', Point2D(textx, texty-.75))
#     print(list(d2d.Offset()))
#     print(list(d2d.GetDrawCoords(Point2D(np.min(crds[:,0]), np.max(crds[:,1])))))
    d2d.FinishDrawing()
    s = d2d.GetDrawingText()
#     s = s.replace('svg:','')
    fig = sg.fromstring(s)
    
    label = sg.TextElement(300, 20, f'{target}: {n1} -> {n2}', size=20, 
                       font='sans-serif', anchor='middle', color='black')
    label2 = sg.TextElement(300, 35, f'{text}', size=15, 
                       font='sans-serif', anchor='middle', color='black')
    label3 = sg.TextElement(150, 285, f'{n1}', size=15, 
                       font='sans-serif', anchor='middle', color='black')
    label4 = sg.TextElement(450, 285, f'{n2}', size=15, 
                       font='sans-serif', anchor='middle', color='black')
    fig.append(label) 
    fig.append(label2) 
    fig.append(label3)
    fig.append(label4)
    return fig.to_str().decode("utf-8") 



def hist(results, c1, c2, ax_max=None, title=''):
    fig = go.Figure()
    nan = [np.isnan(row[c1]) for i, row in results.iterrows()]
    results = results.loc[np.invert(nan), :]
    nan = [np.isnan(row[c2]) for i, row in results.iterrows()]
    results = results.loc[np.invert(nan), :]
    diff = results[c1]-results[c2]

    mu, std = norm.fit(np.array(diff, dtype=float))

    if ax_max is None:
        ax_max=np.max(norm.pdf(np.linspace(-7,7), mu, std)*diff.shape[0])
    fig.update_layout(
        width=600,
        title=title,
        barmode='stack', 
        xaxis = dict(title='DDG(parsley)-DDG(exp) [kcal mol<sup>-1</sup>]', range=(-8,8)),
        yaxis = dict(title='Count', range=(0,ax_max)),
        colorway=colors.qualitative.Safe + colors.qualitative.Vivid,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.01
        )
    )
    max_color = len(fig.layout.colorway)
    for i, target in enumerate(targets.target_dict):
        fig.add_trace(
            go.Histogram(
                x=diff[results['target']==target], 
                xbins_size=0.5, 
                name=target,
                marker_color=fig.layout.colorway[i%max_color],
                legendgroup='group1'
            )
        )
        
    fig.add_trace(
        go.Scatter(
        x=np.linspace(-7,7),
        y=norm.pdf(np.linspace(-7,7), mu, std)*diff.shape[0]*.5,
            line_color='black',
            line_width=3,
            name=f'$\mu={mu:.1f}, \sigma={std:.1f}$'    
        )
    )  
    fig.update_layout(
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.01
        )
    )
    return fig


def create_perturbation_visualization(df, text='', img_size=('400px', '200px'), directory='13_outliers', redraw=False):
    import benchmarkpl
    path = benchmarkpl.__path__[0]
    targets.set_data_dir(path)
    # check whether image exists
    os.makedirs(os.path.join(path, targets.get_target_dir(df["target"]), directory), exist_ok=True)
    file_path = os.path.join(path, targets.get_target_dir(df["target"]), directory, f'{df["edge"]}.svg')
    if not redraw and os.path.exists(file_path):
        with open(file_path, 'r') as file:
            img = file.read()
    else:
        # visualization
        if os.path.exists('../../../02_benchmark_calculations/'):
            targets.set_data_dir('../../../02_benchmark_calculations/')
        target_path = f'{targets.data_path}/{targets.get_target_dir(df["target"])}'
        m1 = Chem.SDMolSupplier(
            f'{target_path}/02_ligands/lig_{df["ligandA"]}/crd/lig_{df["ligandA"]}.sdf', 
            removeHs=False)[0]
        m2 = Chem.SDMolSupplier(
            f'{target_path}/02_ligands/lig_{df["ligandB"]}/crd/lig_{df["ligandB"]}.sdf', 
            removeHs=False)[0]
        pairs = np.loadtxt(
            f'{target_path}/03_hybrid/edge_{df["ligandA"]}_{df["ligandB"]}/water/crd/pairs.dat'
        )
        # decrement pairs to match rdkit counting from 0!
        pairs -= 1
        
        img = drawPerturbationInverted(m1, # rdkit molecule 1
                                       m2, # rdkit molecule 2
                                       pairs, # pairs, np array or list of lists
                                       target=df["target"], # string with target name
                                       n1=df["ligandA"], # name mol 1
                                       n2=df["ligandB"], # name  mol 2
                                       text=text # additional text
                                      )
        
        with open(file_path, 'w') as file:
            file.write(img)
        
        targets.set_data_dir(path)
    original = sg.fromstring(img)
    original.set_size(img_size)
    svgstring = original.to_str().decode("utf-8").rstrip()
    svgstring = '\n'.join(svgstring.split('\n')[1:])
    return svgstring
