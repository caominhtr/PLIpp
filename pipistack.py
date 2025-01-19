import pandas as pd
import numpy as np
from numpy import dot
from numpy.linalg import norm
from scipy.spatial.distance import cdist
import random
import math
from itertools import combinations
import warnings
import sys
import subprocess
import extractfile

file_path = sys.argv[1]

df1, df3, df5, data_list, list_new_indices, list_new = extractfile.process_file(file_path)


def cosine(a, b):
    cos = dot(a,b)/(norm(a)*norm(b))
    return cos

def cosine_pairwise(matrix):
    num_rows = matrix.shape[0]
    similarities = np.zeros((num_rows, num_rows))

    for i in range(num_rows):
        for j in range(num_rows):
            similarities[i, j] = cosine(matrix[i], matrix[j])

    return similarities

def coplanar(mat):
    coplanar_mat = []
    for i in range(0, mat.shape[0]-2):
        coplanar_mat.append(np.dot(np.cross(mat[i], mat[i+1]), mat[i+2])) 
    
    if np.all(np.abs(coplanar_mat) < 0.05):
            return True
    else:
        return False
    
def aromatic_ligand(ligand):

    List = ['C','N','O']
    Result = []
    Center = []
    Mat = []
    Name_ligand = []
    Chain = []

    for res in ligand['Res'].unique():
        for chain in ligand['Chain'].unique():
            selected_data = ligand[(ligand['Name'].isin(List)) & (ligand['Res'] == res) & (ligand['Chain'] == chain)]

            end = len(selected_data)
            n_values = 6
            selected_indices = set()

            data_values = selected_data[['X', 'Y', 'Z']].values

            combs = combinations(range(end), n_values)

            for comb in combs:
                sorted_comb = tuple(sorted(comb))
                selected_indices.add(sorted_comb)

            unique_list = list(selected_indices)

            for uni in np.array(unique_list):
                arrays = data_values[uni]
                center = sum(arrays) / len(arrays)
                vector = (center-arrays).astype('float')
                round_values = np.abs((cosine_pairwise(vector)))
                dist = np.round(np.linalg.norm(vector, axis = 1),2)

                if np.all(np.logical_or(abs(round_values - 0.5) < 10e-2, abs(round_values - 1) < 10e-2)) and np.isclose(dist[0], dist[1:], atol=10e-1).all():
                    if len(selected_data.iloc[uni].Desc.unique()) == len(selected_data.iloc[uni].Desc):
                        if coplanar(vector):
                                Result.append(uni)
                                Center.append(center)
                                Mat.append(data_values[np.array(uni)])
                                Name_ligand.append(res)
                                Chain.append(chain)
    Mat_mat = pd.DataFrame()
    for i in range(len(Mat)):
         Mat_mat = pd.concat([Mat_mat, pd.DataFrame(Mat[i])])
    list_desc = []
    for i in range(int(Mat_mat.shape[0]/6)):
        list_desc.append(list(ligand[(ligand['X'].isin(Mat_mat.iloc[i*6:i*6+6][0])) & (ligand['Y'].isin(Mat_mat.iloc[i*6:i*6+6][1])) & (ligand['Z'].isin(Mat_mat.iloc[i*6:i*6+6][2]))].Desc))

    return Mat_mat, Center, Name_ligand, Chain, list_desc

def aromatic_protein(protein):
    list_PHE = ['CG','CD1','CD2','CE1','CE2','CZ']
    list_TYR = ['CG','CD1','CD2','CE1','CE2','CZ']
    list_TRP = ['CD2','CE2','CE3','CZ2','CZ3','CH2']

    PHE = protein[np.logical_or(np.logical_or(protein['Res']=='APHE', protein['Res'] == 'BPHE'),protein['Res'] == 'PHE') & (protein['Desc'].isin(list_PHE))]
    TYR = protein[np.logical_or(np.logical_or(protein['Res']=='ATYR', protein['Res'] == 'BTYR'),protein['Res'] == 'TYR') & (protein['Desc'].isin(list_TYR))]
    TRP = protein[np.logical_or(np.logical_or(protein['Res']=='ATRP', protein['Res'] == 'BTRP'),protein['Res'] == 'TRP') & (protein['Desc'].isin(list_TRP))]

    aromatic = pd.concat([PHE,TYR,TRP])

    return aromatic

def parallel (matrix_1, matrix_2):
    para = np.array([])
    for i in range(4):
        for j in range(4):
            vect1 = np.cross(matrix_1[i] - matrix_1[i+1], matrix_1[i] - matrix_1[i+2])
            vect2 = np.cross(matrix_2[j] - matrix_2[j+1], matrix_2[j] - matrix_2[j+2])
            para = np.append(para, abs(cosine(vect1, vect2)))

    A = np.sum(np.isclose(para, 1, atol = 0.1))/len(para)
    
    if A > 0.9:
        return True
    else:
        return False 

def perpendicular (matrix_1, matrix_2):
    perpen = np.array([])
    for i in range(4):
        for j in range(4):
            vect1 = np.cross(matrix_1[i] - matrix_1[i+1], matrix_1[i] - matrix_1[i+2])
            vect2 = np.cross(matrix_2[j] - matrix_2[j+1], matrix_2[j] - matrix_2[j+2])
            perpen = np.append(perpen, abs(cosine(vect1, vect2)))

    B = np.sum(np.isclose(perpen, 0, atol = 0.35)) / len(perpen)

    if B > 0.9:
        return True
    else:
        return False
        
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
                            C = ligand.iloc[i][['X','Y','Z']].values
                            O1 = ligand.iloc[j][['X','Y','Z']].values
                            O2 = ligand.iloc[k][['X','Y','Z']].values

                            if np.isclose(cosine(C-O1, C-O2), -0.5, atol = 0.1) and np.isclose(np.linalg.norm(C-O1), 1.2, 0.05) and np.isclose(np.linalg.norm(C-O2), 1.2, 0.05):
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
        if arr_ammonium[0] == 'N' and np.logical_or(len(arr_ammonium) == 5, len(arr_ammonium) == 4):
            list_ammonium.append(data_list[list_new_indices[i_ammonium]][0])

    for i_guanidine, arr_guanidine in enumerate(list_new):
        if np.array_equal(arr_guanidine, [['C','N','N','N']]):
            list_guanidine.append(data_list[list_new_indices[i_guanidine]][0])

    positive_indice = np.concatenate((list_guanidine, list_ammonium, list_sulfonium)).astype('int')
    negative_indice = np.concatenate((list_carboxylate, list_phosphate, list_sulfate)).astype('int')

    positive = ligand[ligand['Order'].isin(positive_indice)][['Order','Chain','Name','Res','X','Y','Z']]
    negative = ligand[ligand['Order'].isin(negative_indice)][['Order','Chain','Name','Res','X','Y','Z']]

    return positive, negative

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
    
def pi_pi_sandwich_stack(matrix_1, matrix_2):

    center_1 = sum(matrix_1) / len(matrix_1)
    center_2 = sum(matrix_2) / len(matrix_2)
    normal_vector = center_1 - center_2
   
    if parallel(matrix_1, matrix_2) and np.linalg.norm(normal_vector) < 5.5:
        return True
    else:
        return False
            
def pi_pi_Tshaped_stack(matrix_1, matrix_2):
    center_1 = sum(matrix_1) / len(matrix_1)
    center_2 = sum(matrix_2) / len(matrix_2)
    normal_vector = center_1 - center_2 

    if perpendicular(matrix_1, matrix_2) and np.linalg.norm(normal_vector) < 5.5:
            return True
    else:
        return False

def pi_pi_interaction(protein, ligand):

    protein_ = aromatic_protein(protein)
    ligand_ = aromatic_ligand(ligand)

    A = ligand_[0].shape[0]
    
    Residue = []
    Number = []
    Chain = []
    Type = []

    Lig_res = []
    Lig_desc = []
    Lig_chain = []
    Lig_type = []

    distance = []

    for index_1 in range(0, protein_.shape[0],6):
        for index_2 in range(0, ligand_[0].shape[0], 6):
            matrix_1 = protein_.iloc[index_1:(index_1 + 6)][['X','Y','Z']].values
            matrix_2 = ligand_[0].iloc[index_2:(index_2 + 6)].values

            center_1 = sum(matrix_1) / len(matrix_1)
            center_2 = sum(matrix_2) / len(matrix_2)
            vector_1 = center_1 - matrix_1
            vector_2 = center_2 - matrix_2

            normal_vector = center_1 - center_2 

            index_2_ = int(index_2 / 6)

            if pi_pi_sandwich_stack(matrix_1, matrix_2):
                Residue.append(protein_.iloc[index_1].Res)
                Number.append(protein_.iloc[index_1].Num)
                Chain.append(protein_.iloc[index_1].Chain)
                Type.append('Pi pi sandwich stack')

                Lig_res.append(ligand_[2][index_2_])
                Lig_desc.append(ligand_[4][index_2_])
                Lig_chain.append(ligand_[3][index_2_])
                Lig_type.append('Pi pi sandwich stack')

                distance.append(np.linalg.norm(normal_vector))



                # print(f'{protein_.iloc[index_1].Res}{protein_.iloc[index_1].Num} chain {protein_.iloc[index_1].Chain} and {ligand_[2][index_2_]} chain {ligand_[3][index_2_]} form a pi_pi sandwich stack with a distance of {round(np.linalg.norm(normal_vector),2)} A')
            if pi_pi_Tshaped_stack(matrix_1, matrix_2):
                Residue.append(protein_.iloc[index_1].Res)
                Number.append(protein_.iloc[index_1].Num)
                Chain.append(protein_.iloc[index_1].Chain)
                Type.append('Pi pi Tshaped stack')


                Lig_res.append(ligand_[2][index_2_])
                Lig_desc.append(ligand_[4][index_2_])
                Lig_chain.append(ligand_[3][index_2_])
                Lig_type.append('Pi pi Tshaped stack')
                
                distance.append(np.linalg.norm(normal_vector))

                # print(f'{protein_.iloc[index_1].Res}{protein_.iloc[index_1].Num} chain {protein_.iloc[index_1].Chain} and {ligand_[2][index_2_]} chain {ligand_[3][index_2_]} form a T-shaped stack with a distance of {round(np.linalg.norm(normal_vector),2)} A')
    
    aro_ligand_ =  ligand_
    aro_ligand = aro_ligand_[0]
    cation_protein = charge_protein(protein)[0]
    
    aro_protein = aromatic_protein(protein)
    cation_ligand = charge_ligand(ligand)[0]
    

    for i in range(int(aro_ligand.shape[0]/6)):
        for j in range(cation_protein.shape[0]):
            cation_pi_angle = []

            aro = aro_ligand.iloc[(i*6):(i*6+6)]
            cation = cation_protein.iloc[j][['X','Y','Z']].values
            center = aro.sum()/6

            R = norm(cation - center)
            
            for k in range(5):
                vect = aro.iloc[k] - aro.iloc[k+1]
                cation_pi_angle.append(cosine(cation-center, vect))
            
            num =  np.sum((np.array(cation_pi_angle) < 0.707) & (np.array(cation_pi_angle) > -0.707))

            if R < 5.5 and num >= 4:
                aro_ligand_res = ligand[ligand['X'] == aro_ligand.iloc[i*6].iloc[0]].Res.values[0]
                aro_ligand_desc = ligand[ligand['X'] == aro_ligand.iloc[i*6].iloc[0]].Desc.values[0]
                aro_ligand_chain = ligand[ligand['X'] == aro_ligand.iloc[i*6].iloc[0]].Chain.values[0]

                cation_protein_res = protein[protein['X'] == cation_protein.iloc[j]['X']].Res.values[0]
                cation_protein_num = protein[protein['X'] == cation_protein.iloc[j]['X']].Num.values[0]
                cation_protein_chain = protein[protein['X'] == cation_protein.iloc[j]['X']].Chain.values[0]

                Residue.append(cation_protein_res)
                Number.append(cation_protein_num)
                Chain.append(cation_protein_chain)
                Type.append('Cation-pi Cation protein')

                Lig_res.append(ligand[ligand['X'] == aro_ligand.iloc[i*6].iloc[0]].Res.values[0])
                Lig_desc.append(aro_ligand_[4][i])
                Lig_chain.append(ligand[ligand['X'] == aro_ligand.iloc[i*6].iloc[0]].Chain.values[0])
                Lig_type.append('Cation-pi aromatic ligand')

                distance.append(R)

                # print(f'Protein (cation): {cation_protein_res} {cation_protein_num}{cation_protein_chain}, Ligand (aromatic): {aro_ligand_res} {aro_ligand_desc} in chain {aro_ligand_chain}')

    for i in range(int(aro_protein.shape[0]/6)):
        for j in range(cation_ligand.shape[0]):
            cation_pi_angle = []

            aro = aro_protein.iloc[(i*6):(i*6+6)][['X','Y','Z']]
            cation = cation_ligand.iloc[j][['X','Y','Z']]
            center = aro.sum()/6

            R = norm(cation - center)
            
            for k in range(5):
                vect = aro.iloc[k] - aro.iloc[k+1]
                cation_pi_angle.append(cosine(cation-center, vect))
            
            num =  np.sum((np.array(cation_pi_angle) < 0.707) & (np.array(cation_pi_angle) > -0.707))

            if R < 5.5 and num >= 4:
                aro_protein_res = protein[protein['X'] == aro.iloc[0].iloc[0]].Res.values[0]
                aro_protein_num = protein[protein['X'] == aro.iloc[0].iloc[0]].Num.values[0]
                aro_protein_chain = protein[protein['X'] == aro.iloc[0].iloc[0]].Chain.values[0]

                cation_ligand_res = ligand[ligand['X'] == cation_ligand.iloc[j]['X']].Res.values[0]
                cation_ligand_desc = ligand[ligand['X'] == cation_ligand.iloc[j]['X']].Desc.values[0]
                cation_ligand_chain = ligand[ligand['X'] == cation_ligand.iloc[j]['X']].Chain.values[0]

                Residue.append(aro_protein_res)
                Number.append(aro_protein_num)
                Chain.append(aro_protein_chain)
                Type.append('Cation-pi Aromatic protein')

                Lig_res.append(ligand[ligand['X'] == cation_ligand.iloc[j]['X']].Res.values[0])
                Lig_desc.append(ligand[ligand['X'] == cation_ligand.iloc[j]['X']].Desc.values[0])
                Lig_chain.append(ligand[ligand['X'] == cation_ligand.iloc[j]['X']].Chain.values[0])
                Lig_type.append('Cation-pi cation ligand')

                distance.append(R)
                # print(f'Protein (aromatic): {aro_protein_res} {aro_protein_num}{aro_protein_chain}, Ligand (cation): {cation_ligand_res} {cation_ligand_desc} in chain {cation_ligand_chain}')


    columns = ['Residue', 'Number', 'Chain','Type']
    df_result = pd.DataFrame(columns=columns)
    df_result['Residue'] = Residue
    df_result['Number'] = Number
    df_result['Chain'] = Chain
    df_result['Type'] = Type
    df_result['Ligand'] = Lig_res
    df_result['Atom'] = Lig_desc
    df_result['Distance'] = distance

    columns_ = ['Ligand', 'Atom', 'Chain','Type']
    df_pharmacophore = pd.DataFrame(columns=columns_)
    df_pharmacophore['Ligand'] = Lig_res
    df_pharmacophore['Atom'] = Lig_desc
    df_pharmacophore['Chain'] = Lig_chain
    df_pharmacophore['Type'] = Lig_type
    df_pharmacophore = df_pharmacophore.drop_duplicates(subset='Atom')
    df_pharmacophore = df_pharmacophore.reset_index(drop=True, inplace=False)

    return (df_result, df_pharmacophore)

# warnings.filterwarnings('ignore', 'Boolean Series key will be reindexed to match DataFrame index')
# print(pi_pi_interaction(df3, df1))

