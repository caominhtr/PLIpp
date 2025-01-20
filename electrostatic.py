import pandas as pd
import numpy as np
from numpy import dot
from numpy.linalg import norm
import warnings
import sys
import extractfile

file_path = sys.argv[1]

df1, df3, df5, data_list, list_new_indices, list_new = extractfile.process_file(file_path)

def cosine(a, b):
    cos = dot(a,b)/(norm(a)*norm(b))
    return cos


def charge_protein (protein):
    asp = protein[np.logical_or(np.logical_or(protein['Res']=='AASP', protein['Res'] == 'BASP'),protein['Res'] == 'ASP') & np.logical_or(np.logical_or(protein['Desc'] == 'CG', protein['Desc'] == 'OD1'), protein['Desc'] == 'OD2')]
    glu = protein[np.logical_or(np.logical_or(protein['Res']=='AGLU', protein['Res'] == 'BGLU'),protein['Res'] == 'GLU') & np.logical_or(np.logical_or(protein['Desc'] == 'CD', protein['Desc'] == 'OE1'), protein['Desc'] == 'OE2')]
    lys = protein[np.logical_or(np.logical_or(protein['Res']=='ALYS', protein['Res'] == 'BLYS'),protein['Res'] == 'LYS') & (protein['Desc'] == 'NZ')]
    arg = protein[np.logical_or(np.logical_or(protein['Res']=='AARG', protein['Res'] == 'BARG'),protein['Res'] == 'ARG') & np.logical_or(np.logical_or(protein['Desc'] == 'CZ', protein['Desc'] == 'NH1'), protein['Desc'] == 'NH2')]
    his = protein[np.logical_or(np.logical_or(protein['Res']=='AHIS', protein['Res'] == 'BHIS'),protein['Res'] == 'HIS') & (protein['Desc'] == 'NE2')]

    asp_center = asp.groupby(['Num','Chain','Res'])[['X','Y','Z']].mean().reset_index()
    glu_center = glu.groupby(['Num','Chain','Res'])[['X','Y','Z']].mean().reset_index()
    arg_center = arg.groupby(['Num','Chain','Res'])[['X','Y','Z']].mean().reset_index()
    lys_ = lys[['Num','Chain','Res','X','Y','Z']].reset_index(drop = True)
    his_ = his[['Num','Chain','Res','X','Y','Z']].reset_index(drop = True)

    positive = pd.concat([lys_, arg_center, his_]).reset_index(drop = True)
    negative = pd.concat([asp_center, glu_center]).reset_index(drop = True)

    return positive, negative

def charge_ligand(ligand):    

    list_carboxylate = []
    list_phosphate = []
    list_sulfate = []
    list_sulfonium = []
    list_ammonium = []
    list_guanidine = []

    for res in ligand['Res'].unique():
        for chain in ligand['Chain'].unique():
            ligand_ = ligand[(ligand['Res'] == res) & (ligand['Chain'] == chain)]
            for i in range(ligand_.shape[0]):
                for j in range(ligand_.shape[0]):
                    for k in range(ligand_.shape[0]):
                        if (ligand_.iloc[i]['Name'] == 'C') and (ligand_.iloc[j]['Name'] == 'O') and (ligand_.iloc[k]['Name'] == 'O'):
                            C = ligand_.iloc[i][['X','Y','Z']].values
                            O1 = ligand_.iloc[j][['X','Y','Z']].values
                            O2 = ligand_.iloc[k][['X','Y','Z']].values

                            if np.isclose(cosine(C-O1, C-O2), -0.5, atol = 0.1) and np.isclose(np.linalg.norm(C-O1), 1.3, atol=0.05) and np.isclose(np.linalg.norm(C-O2), 1.3, atol=0.05):
                                list_carboxylate.append(ligand_.iloc[i].Order)
    list_carboxylate = list(set(list_carboxylate))


    for i_phosphate, arr_phosphate in enumerate(list_new):
        if arr_phosphate[0] == 'P' and len(arr_phosphate) == 5:
            list_phosphate.append(data_list[list_new_indices[i_phosphate]][0])

    for i_sulfate, arr_sulfate in enumerate(list_new):
        if arr_sulfate[0] == 'S' and len(arr_sulfate) == 5:
            list_sulfate.append(data_list[list_new_indices[i_sulfate]][0])

    for i_sulfonium, arr_sulfonium in enumerate(list_new):
        if arr_sulfonium[0] == 'S' and len(arr_sulfonium) == 4:
            list_sulfonium.append(data_list[list_new_indices[i_sulfonium]][0])

    for i_ammonium, arr_ammonium in enumerate(list_new):
        if arr_ammonium[0] == 'N' and np.logical_or(len(arr_ammonium) == 5, len(arr_ammonium) == 4) and np.all(arr_ammonium != 'H'):
            list_ammonium.append(data_list[list_new_indices[i_ammonium]][0])

    for i_guanidine, arr_guanidine in enumerate(list_new):
        if np.array_equal(arr_guanidine, [['C','N','N','N']]):
            list_guanidine.append(data_list[list_new_indices[i_guanidine]][0])

    positive_indice = np.concatenate((list_guanidine, list_ammonium, list_sulfonium)).astype('int')
    negative_indice = np.concatenate((list_carboxylate, list_phosphate, list_sulfate)).astype('int')

    positive = ligand[ligand['Order'].isin(positive_indice)][['Order','Chain','Name','Res','X','Y','Z', 'Desc']]
    negative = ligand[ligand['Order'].isin(negative_indice)][['Order','Chain','Name','Res','X','Y','Z', 'Desc']]

    return positive, negative

def electrostatics_interaction(protein, ligand):
    pro_charge_positive = charge_protein(protein)[0]
    pro_charge_negative = charge_protein(protein)[1]

    lig_charge_positive = charge_ligand(ligand)[0]
    lig_charge_negative = charge_ligand(ligand)[1]

    Residue = []
    Number = []
    Chain = []
    Type = []

    Lig_res = []
    Lig_desc = []
    Lig_chain = []
    Lig_type = []

    distance = []

    for pro_pos in range(pro_charge_positive.shape[0]):
        for lig_neg in range(lig_charge_negative.shape[0]):
            if norm(pro_charge_positive.iloc[pro_pos][['X','Y','Z']] - lig_charge_negative.iloc[lig_neg][['X','Y','Z']]) < 5.5:
                Residue.append(pro_charge_positive.iloc[pro_pos]['Res'])
                Number.append(pro_charge_positive.iloc[pro_pos]['Num'])
                Chain.append(pro_charge_positive.iloc[pro_pos]['Chain'])
                Type.append('Electrostatic positive')

                distance.append(norm(pro_charge_positive.iloc[pro_pos][['X','Y','Z']] - lig_charge_negative.iloc[lig_neg][['X','Y','Z']]))

                Lig_res.append(lig_charge_negative.iloc[lig_neg]['Res'])
                Lig_desc.append(lig_charge_negative.iloc[lig_neg]['Desc'])
                Lig_chain.append(lig_charge_negative.iloc[lig_neg]['Chain'])
                Lig_type.append('Electrostatic negative')
                # print(f'{pro_charge_positive.iloc[pro_pos]['Res'],pro_charge_positive.iloc[pro_pos]['Num'],pro_charge_positive.iloc[pro_pos]['Chain']} and {lig_charge_negative.iloc[lig_neg]['Name'],lig_charge_negative.iloc[lig_neg]['Res'], lig_charge_negative.iloc[lig_neg]['Chain']} form an electrostatic bond')

    for pro_neg in range(pro_charge_negative.shape[0]):
        for lig_pos in range(lig_charge_positive.shape[0]):
            if norm(pro_charge_negative.iloc[pro_neg][['X','Y','Z']] - lig_charge_positive.iloc[lig_pos][['X','Y','Z']]) < 5.5:
                Residue.append(pro_charge_negative.iloc[pro_neg]['Res'])
                Number.append(pro_charge_negative.iloc[pro_neg]['Num'])
                Chain.append(pro_charge_negative.iloc[pro_neg]['Chain'])
                Type.append('Electrostatic negative')
                distance.append(norm(pro_charge_negative.iloc[pro_neg][['X','Y','Z']] - lig_charge_positive.iloc[lig_pos][['X','Y','Z']]))

                Lig_res.append(lig_charge_positive.iloc[lig_pos]['Res'])
                Lig_desc.append(lig_charge_positive.iloc[lig_pos]['Desc'])
                Lig_chain.append(lig_charge_positive.iloc[lig_pos]['Chain'])
                Lig_type.append('Electrostatic positive')
                # print(f'{pro_charge_negative.iloc[pro_neg]['Res'],pro_charge_negative.iloc[pro_neg]['Num'], pro_charge_positive.iloc[pro_pos]['Chain']} and {lig_charge_positive.iloc[lig_pos]['Name'],lig_charge_positive.iloc[lig_pos]['Res'], lig_charge_positive.iloc[lig_pos]['Chain']} form an electrostatic bond')


    columns = ['Residue', 'Number', 'Chain','Type']
    df_result = pd.DataFrame(columns=columns)
    df_result['Residue'] = Residue
    df_result['Number'] = Number
    df_result['Chain'] = Chain
    df_result['Type'] = Type
    df_result['Ligand'] = Lig_res
    df_result['Atom'] = Lig_desc
    df_result['Distance'] = distance
    df_result = df_result.drop_duplicates()
    df_result = df_result.reset_index(drop=True, inplace=False)

    columns_ = ['Ligand', 'Atom', 'Chain','Type']
    df_pharmacophore = pd.DataFrame(columns=columns_)
    df_pharmacophore['Ligand'] = Lig_res
    df_pharmacophore['Atom'] = Lig_desc
    df_pharmacophore['Chain'] = Lig_chain
    df_pharmacophore['Type'] = Lig_type
    df_pharmacophore = df_pharmacophore.drop_duplicates()
    df_pharmacophore = df_pharmacophore.reset_index(drop=True, inplace=False)

    return df_result, df_pharmacophore

# print(electrostatics_interaction(df3, df1)[0])
# print(electrostatics_interaction(df3, df1)[1])
