import pandas as pd
import numpy as np
from numpy.linalg import norm
import warnings
import sys
import extractfile

file_path = sys.argv[1]

df1, df3, df5, data_list, list_new_indices, list_new = extractfile.process_file(file_path)

def hydrophobic_interaction(protein, ligand):

    list = [('CZ','TYR'), ('CZ', 'ATYR'), ('CZ', 'BTYR'),
            ('CB','SER'), ('CB', 'ASER'), ('CB', 'BSER'),
            ('CB','THR'), ('CB', 'ATHR'), ('CB', 'BTHR'),
            ('CB','CYS'), ('CB', 'ACYS'), ('CB', 'BCYS'),
            ('CG','MET'),('CE','MET'), ('CG', 'AMET'), ('CE', 'AMET'), ('CG', 'BMET'), ('CE', 'BMET'),
            ('CD','GLN'), ('CD', 'AGLN'), ('CD', 'BGLN'),
            ('CG','ASP'), ('CG', 'AASP'), ('CG', 'BASP'),
            ('CD','GLU'), ('CD', 'AGLU'), ('CD', 'BGLU'),
            ('CZ','ARG'), ('CZ', 'AARG'), ('CZ', 'BARG'),
            ('CE','LYS'), ('CE', 'ALYS'), ('CE', 'BLYS'),
            ('CD','PRO'), ('CD', 'APRO'), ('CD', 'BPRO'),
            ('CD1', 'TRP'), ('CE2', 'TRP'), ('CD1', 'ATRP'), ('CE2', 'ATRP'), ('CD1', 'BTRP'), ('CE2', 'BTRP'),
            ('CG','HIS'),('CD2','HIS'),('CE1','HIS'), ('CG', 'AHIS'), ('CD2', 'AHIS'), ('CE1', 'AHIS'), ('CG', 'BHIS'), ('CD2', 'BHIS'), ('CE1', 'BHIS'),
            ('CG','ASN'),('CD','GLN'), ('CG', 'AASN'), ('CD', 'AASN'), ('CG', 'BASN'), ('CD', 'BASN')]
    
    protein_ = protein[protein['Name'] == 'C']
    ligand_ = ligand[ligand['Name'] == 'C'] 

    Carbon_remove = []

    for i_carbon, arr_carbon in enumerate(list_new):
        if arr_carbon[0] == 'N' and arr_carbon[1] == 'C':
            Carbon_remove.append(data_list[list_new_indices[i_carbon]][1])

    for i_carbon, arr_carbon in enumerate(list_new):
        if arr_carbon[0] == 'O' and arr_carbon[1] == 'C':
            Carbon_remove.append(data_list[list_new_indices[i_carbon]][1])

    for i_carbon, arr_carbon in enumerate(list_new):
        if arr_carbon[0] == 'C' and np.logical_or(np.logical_or(arr_carbon[1] == 'N', arr_carbon[1] == 'O'), arr_carbon[1] == 'S'):
            Carbon_remove.append(data_list[list_new_indices[i_carbon]][0])

    for i_carbon, arr_carbon in enumerate(list_new):
        if arr_carbon[0] == 'C' and np.logical_or(np.logical_or(np.logical_or(arr_carbon[1] == 'F', arr_carbon[1] == 'Cl'), arr_carbon[1] == 'Br'), arr_carbon[1] == 'I'):
            Carbon_remove.append(data_list[list_new_indices[i_carbon]][0])

    for i_carbon, arr_carbon in enumerate(list_new):
        if arr_carbon[1] == 'C' and np.logical_or(np.logical_or(np.logical_or(arr_carbon[0] == 'F', arr_carbon[0] == 'Cl'), arr_carbon[0] == 'Br'), arr_carbon[0] == 'I'):
            Carbon_remove.append(data_list[list_new_indices[i_carbon]][1])

    for i_carbon, arr_carbon in enumerate(list_new):
        if arr_carbon[0] == 'C' and np.logical_or(any(arr_carbon == 'N'), any(arr_carbon == 'O')):
            Carbon_remove.append(data_list[list_new_indices[i_carbon]][0])

    Residue = []
    Number = []
    Chain = []

    Lig_res = []
    Lig_desc = []
    Lig_chain = []

    distance = []

    for index_1 in range(protein_.shape[0]):
        for index_2 in range(ligand_.shape[0]):
            if protein_.iloc[index_1].Desc not in ['C','CA']:
                if (protein_.iloc[index_1].Desc, protein_.iloc[index_1].Res) not in list:
                    if ligand_.iloc[index_2].Order not in Carbon_remove:
                        dist_1 = protein_.iloc[index_1][['X','Y','Z']].values
                        dist_2 = ligand_.iloc[index_2][['X','Y','Z']].values
                        dist = np.linalg.norm(dist_1 - dist_2)

                        if dist >= 3.3 and dist <= 4.0:
                            Residue.append(protein_.iloc[index_1].Res)
                            Number.append(protein_.iloc[index_1].Num)
                            Chain.append(protein_.iloc[index_1].Chain)

                            Lig_res.append(ligand_.iloc[index_2].Res)
                            Lig_desc.append(ligand_.iloc[index_2].Desc)
                            Lig_chain.append(ligand_.iloc[index_2].Chain)

                            distance.append(dist)

                            # print(f' {protein_.iloc[index_1].Desc, protein_.iloc[index_1].Res, protein_.iloc[index_1].Chain,protein_.iloc[index_1].Num} and {ligand_.iloc[index_2].Res ,ligand_.iloc[index_2].Desc, ligand_.iloc[index_2].Chain} form a hydrophobic interaction')
    
    columns = ['Residue', 'Number', 'Chain','Type']
    df_result = pd.DataFrame(columns=columns)
    df_result['Residue'] = Residue
    df_result['Number'] = Number
    df_result['Chain'] = Chain
    df_result['Type'] = np.repeat('Hydrophobic', len(df_result))
    df_result['Ligand'] = Lig_res
    df_result['Atom'] = Lig_desc
    df_result['Distance'] = distance
    df_result = df_result.drop_duplicates(keep='first')
    df_result = df_result.reset_index(drop=True, inplace=False)

    columns_ = ['Ligand', 'Atom', 'Chain','Type']
    df_pharmacophore = pd.DataFrame(columns=columns_)
    df_pharmacophore['Ligand'] = Lig_res
    df_pharmacophore['Atom'] = Lig_desc
    df_pharmacophore['Chain'] = Lig_chain
    df_pharmacophore['Type'] = np.repeat('Hydrophobic', len(df_pharmacophore))
    df_pharmacophore = df_pharmacophore.drop_duplicates(keep='first')
    df_pharmacophore = df_pharmacophore.reset_index(drop=True, inplace=False)

    return df_result, df_pharmacophore

# print(hydrophobic_interaction(df3, df1)[0])
# print(hydrophobic_interaction(df3, df1)[1])
