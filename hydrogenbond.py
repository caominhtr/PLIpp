import pandas as pd
import numpy as np
from numpy import dot
from numpy.linalg import norm
from scipy.spatial.distance import cdist
import warnings
import sys
import extractfile

file_path = sys.argv[1]

df1, df3, df5, data_list, list_new_indices, list_new = extractfile.process_file(file_path)

def hydrogen_bond_protein(protein):
    ser = protein[np.logical_or(np.logical_or(protein['Res']=='ASER', protein['Res'] == 'BSER'),protein['Res'] == 'SER') & (protein['Desc'] == 'OG')]
    thr = protein[np.logical_or(np.logical_or(protein['Res']=='ATHR', protein['Res'] == 'BTHR'),protein['Res'] == 'THR') & (protein['Desc'] == 'OG1')]
    tyr = protein[np.logical_or(np.logical_or(protein['Res']=='ATYR', protein['Res'] == 'BTYR'),protein['Res'] == 'TYR') & (protein['Desc'] == 'OH')]
    arg = protein[np.logical_or(np.logical_or(protein['Res']=='AARG', protein['Res'] == 'BARG'),protein['Res'] == 'ARG') & (np.logical_or(protein['Desc'] == 'NH1', protein['Desc'] == 'NH2'))]
    lys = protein[np.logical_or(np.logical_or(protein['Res']=='ALYS', protein['Res'] == 'BLYS'),protein['Res'] == 'LYS') & (protein['Desc'] == 'NZ')]
    gln = protein[np.logical_or(np.logical_or(protein['Res']=='AGLN', protein['Res'] == 'BGLN'),protein['Res'] == 'GLN') & (protein['Desc'] == 'NE2')]
    asn = protein[np.logical_or(np.logical_or(protein['Res']=='AASN', protein['Res'] == 'BASN'),protein['Res'] == 'ASN') & (protein['Desc'] == 'ND2')]
    his = protein[np.logical_or(np.logical_or(protein['Res']=='AHIS', protein['Res'] == 'BHIS'),protein['Res'] == 'HIS') & (protein['Desc'] == 'ND1')]
    trp = protein[np.logical_or(np.logical_or(protein['Res']=='ATRP', protein['Res'] == 'BTRP'),protein['Res'] == 'TRP') & (protein['Desc'] == 'NE1')]
    cys = protein[np.logical_or(np.logical_or(protein['Res']=='ACYS', protein['Res'] == 'BCYS'),protein['Res'] == 'CYS') & (protein['Desc'] == 'SG')]
    nitrogen = protein[(protein['Desc'] == 'N')]

    ser_H = protein[np.logical_or(np.logical_or(protein['Res']=='ASER', protein['Res'] == 'BSER'),protein['Res'] == 'SER') & (protein['Desc'] == 'HG')]
    thr_H = protein[np.logical_or(np.logical_or(protein['Res']=='ATHR', protein['Res'] == 'BTHR'),protein['Res'] == 'THR') & (protein['Desc'] == 'HG1')]
    tyr_H = protein[np.logical_or(np.logical_or(protein['Res']=='ATYR', protein['Res'] == 'BTYR'),protein['Res'] == 'TYR') & (protein['Desc'] == 'HH')]
    arg_H1 = protein[np.logical_or(np.logical_or(protein['Res']=='AARG', protein['Res'] == 'BARG'),protein['Res'] == 'ARG') & (np.logical_or(protein['Desc'] == 'HH11', protein['Desc'] == 'HH12'))]
    arg_H2 = protein[np.logical_or(np.logical_or(protein['Res']=='AARG', protein['Res'] == 'BARG'),protein['Res'] == 'ARG') & (np.logical_or(protein['Desc'] == 'HH21', protein['Desc'] == 'HH22'))]
    lys_H = protein[np.logical_or(np.logical_or(protein['Res']=='ALYS', protein['Res'] == 'BLYS'),protein['Res'] == 'LYS') & (np.logical_or(np.logical_or(protein['Desc'] == 'HZ1', protein['Desc'] == 'HZ2'), protein['Desc'] == 'HZ3'))]
    gln_H = protein[np.logical_or(np.logical_or(protein['Res']=='AGLN', protein['Res'] == 'BGLN'),protein['Res'] == 'GLN') & (np.logical_or(protein['Desc'] == 'HE21', protein['Desc'] == 'HE22'))]
    asn_H = protein[np.logical_or(np.logical_or(protein['Res']=='AASN', protein['Res'] == 'BASN'),protein['Res'] == 'ASN') & (np.logical_or(protein['Desc'] == 'HD21', protein['Desc'] == 'HD22'))]
    his_H = protein[np.logical_or(np.logical_or(protein['Res']=='AHIS', protein['Res'] == 'BHIS'),protein['Res'] == 'HIS') & (protein['Desc'] == 'HD1')]
    trp_H = protein[np.logical_or(np.logical_or(protein['Res']=='ATRP', protein['Res'] == 'BTRP'),protein['Res'] == 'TRP') & (protein['Desc'] == 'HE1')]
    cys_H = protein[np.logical_or(np.logical_or(protein['Res']=='ACYS', protein['Res'] == 'BCYS'),protein['Res'] == 'CYS') & (protein['Desc'] == 'HG')]
    nitrogen_H = protein[(protein['Desc'] == 'H')]

    carbonyl = protein[protein['Desc'] == 'O']
    met = protein[np.logical_or(np.logical_or(protein['Res']=='AMET', protein['Res'] == 'BMET'),protein['Res'] == 'MET') & (protein['Desc'] == 'SD')]
    carbonyl_asp = protein[np.logical_or(np.logical_or(protein['Res']=='AASP', protein['Res'] == 'BASP'),protein['Res'] == 'ASP') & (np.logical_or(protein['Desc'] == 'OD1', protein['Desc'] == 'OD2'))]
    carbonyl_glu = protein[np.logical_or(np.logical_or(protein['Res']=='AGLU', protein['Res'] == 'BGLU'),protein['Res'] == 'GLU') & (np.logical_or(protein['Desc'] == 'OE1', protein['Desc'] == 'OE2'))]
    carbonyl_asn = protein[np.logical_or(np.logical_or(protein['Res']=='AASN', protein['Res'] == 'BASN'),protein['Res'] == 'ASN') & (protein['Desc'] == 'OD1')]
    carbonyl_gln = protein[np.logical_or(np.logical_or(protein['Res']=='AGLN', protein['Res'] == 'BGLN'),protein['Res'] == 'GLN') & (protein['Desc'] == 'OE1')]

    donor = pd.concat([ser, thr, tyr,arg, lys, gln, asn, his, trp, cys, nitrogen])
    hydrogen = pd.concat([ser_H, thr_H, tyr_H,  arg_H1, arg_H2 ,lys_H, gln_H, asn_H, his_H, trp_H, cys_H ,nitrogen_H])
    acceptor = pd.concat([donor, carbonyl, met, carbonyl_asp, carbonyl_glu, carbonyl_asn, carbonyl_gln])


    return donor, hydrogen, acceptor

def hydrogen_bond_ligand(ligand):

    oxygen = ligand[ligand['Name'] == 'O']
    nitrogen = ligand[ligand['Name'] == 'N']
    sulfur = ligand[ligand['Name'] == 'S']

    hydrogen = ligand[ligand['Name'] == 'H']

    chloro = ligand[ligand['Name'] == 'Cl']
    bromo = ligand[ligand['Name'] == 'Br']
    fluoro = ligand[ligand['Name'] == 'F']
    iodo = ligand[ligand['Name'] == 'I']

    donor = pd.concat([oxygen, nitrogen, sulfur])
    acceptor = pd.concat([donor, chloro, bromo, fluoro, iodo])

    return donor, hydrogen, acceptor

def cosine(a, b):
    cos = dot(a,b)/(norm(a)*norm(b))
    return cos

def hydrogen_bond_interaction(protein, ligand):
    data_protein = hydrogen_bond_protein(protein)
    data_ligand = hydrogen_bond_ligand(ligand)

    distances_protein_donor_hydrogen = cdist(data_protein[0][['X','Y','Z']], data_protein[1][['X','Y','Z']])
    distances_ligand_donor_hydrogen = cdist(data_ligand[0][['X','Y','Z']], data_ligand[1][['X','Y','Z']])


    indices_1= np.where(distances_protein_donor_hydrogen < 1.1)
    indices_2 = np.where(distances_ligand_donor_hydrogen < 1.1)


    protein_donor = data_protein[0].iloc[indices_1[0]]
    protein_hydrogen = data_protein[1].iloc[indices_1[1]]
    protein_acceptor = data_protein[2]

    ligand_donor = data_ligand[0].iloc[indices_2[0]]
    ligand_hydrogen = data_ligand[1].iloc[indices_2[1]]
    ligand_acceptor = data_ligand[2]

    dist_1 = cdist(protein_donor[['X','Y','Z']], ligand_acceptor[['X','Y','Z']])
    dist_2 = cdist(protein_acceptor[['X','Y','Z']], ligand_donor[['X','Y','Z']])

    indices_3 = np.where(np.logical_and((dist_1 > 2.5),(dist_1 < 3.8)))
    indices_4 = np.where(np.logical_and((dist_2 > 2.5),(dist_2 < 3.8)))

    protein_donor_ = protein_donor.iloc[indices_3[0]].drop_duplicates()
    protein_hydrogen_ = protein_hydrogen.iloc[indices_3[0]].drop_duplicates()
    protein_acceptor_ = protein_acceptor.iloc[indices_4[0]].drop_duplicates()

    ligand_donor_ = ligand_donor.iloc[indices_4[1]].drop_duplicates()
    ligand_hydrogen_ = ligand_hydrogen.iloc[indices_4[1]].drop_duplicates()
    ligand_acceptor_ = ligand_acceptor.iloc[indices_3[1]].drop_duplicates()

    Residue = []
    Number = []
    Chain = []
    Type = []

    Lig_res = []
    Lig_desc = []
    Lig_chain = []
    Lig_type = []

    distance = []


    for donor_1 in range(protein_donor_.shape[0]):
        for hydro_1 in range(protein_hydrogen_.shape[0]):
            for acceptor_1 in range(ligand_acceptor_.shape[0]):

                D_1 = protein_donor_.iloc[donor_1][['X','Y','Z']]
                H_1 = protein_hydrogen_.iloc[hydro_1][['X','Y','Z']]
                A_1 = ligand_acceptor_.iloc[acceptor_1][['X','Y','Z']]

                if cosine(H_1-D_1, H_1-A_1) < -0.65 and norm(D_1 - H_1) < 1.1 and norm(D_1 - A_1) > 2.5 and norm(D_1 - A_1) < 3.8 :
                    Residue.append(protein_donor_.iloc[donor_1].Res)
                    Number.append(protein_donor_.iloc[donor_1].Num)
                    Chain.append(protein_donor_.iloc[donor_1].Chain)
                    Type.append('Hydrogen bond donor')

                    Lig_res.append(ligand_acceptor_.iloc[acceptor_1].Res)
                    Lig_desc.append(ligand_acceptor_.iloc[acceptor_1].Desc)
                    Lig_chain.append(ligand_acceptor_.iloc[acceptor_1].Chain)
                    Lig_type.append('Hydrogen bond acceptor')

                    distance.append(norm(D_1 - A_1))
                    # print(f'{protein_donor_.iloc[donor_1].Res, protein_donor_.iloc[donor_1].Num, protein_donor_.iloc[donor_1].Desc, protein_hydrogen_.iloc[hydro_1].Desc} chain {protein_donor_.iloc[donor_1].Chain} and {ligand_acceptor_.iloc[acceptor_1].Res, ligand_acceptor_.iloc[acceptor_1].Desc} form a hydrogen bond with {cosine(H_1-D_1, H_1-A_1)}, {norm(H_1-D_1)}, {norm(H_1-A_1)}, {norm(D_1-A_1)}')

    
    for donor_2 in range(ligand_donor_.shape[0]):
        for hydro_2 in range(ligand_hydrogen_.shape[0]):
            for acceptor_2 in range(protein_acceptor_.shape[0]):

                D_2 = ligand_donor_.iloc[donor_2][['X','Y','Z']]
                H_2 = ligand_hydrogen_.iloc[hydro_2][['X','Y','Z']]
                A_2 = protein_acceptor_.iloc[acceptor_2][['X','Y','Z']]

                if cosine(H_2-D_2, H_2-A_2) < -0.65 and norm(D_2 - H_2) < 1.1 and norm(D_2 - A_2) > 2.5 and norm(D_2 - A_2) < 3.8:
                    Residue.append(protein_acceptor_.iloc[acceptor_2].Res)
                    Number.append(protein_acceptor_.iloc[acceptor_2].Num)
                    Chain.append(protein_acceptor_.iloc[acceptor_2].Chain)
                    Type.append('Hydrogen bond acceptor')

                    Lig_res.append(ligand_donor_.iloc[donor_2].Res)
                    Lig_desc.append(ligand_donor_.iloc[donor_2].Desc)
                    Lig_chain.append(ligand_donor_.iloc[donor_2].Chain)
                    Lig_type.append('Hydrogen bond donor')

                    distance.append(norm(D_2 - A_2))
                    # print(f'{ligand_donor_.iloc[donor_2].Res, ligand_donor_.iloc[donor_2].Desc, ligand_hydrogen_.iloc[hydro_2].Desc} and {protein_acceptor_.iloc[acceptor_2].Res, protein_acceptor_.iloc[acceptor_2].Num, protein_acceptor_.iloc[acceptor_2].Desc} chain {protein_acceptor_.iloc[acceptor_2].Chain} form a hydrogen bond with {cosine(H_2-D_2, H_2-A_2)}, {norm(H_2-D_2)}, {norm(H_2-A_2)}, {norm(D_2-A_2)}')
    

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

# print(hydrogen_bond_interaction(df3, df1)[0], hydrogen_bond_interaction(df3, df1)[1])