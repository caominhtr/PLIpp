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

def halogen_bond_interaction(protein, ligand):
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
    halogen_list = ["Cl", "Br", "I"]

    halogen = ligand[ligand["Name"].isin(halogen_list)][['X','Y','Z']]
    halogen_donor = ligand[ligand["Name"] == "C"][['X','Y','Z']]
    halogen_acceptor = hydrogen_bond_protein(protein)[2][['X','Y','Z']]

    Residue = []
    Number = []
    Chain = []
    Type = []
    distance = []

    Lig_res = []
    Lig_num = []
    Lig_chain = []

    for i in range(halogen.shape[0]):
        for j in range(halogen_donor.shape[0]):
            for k in range(halogen_acceptor.shape[0]):
                length_bond_halogen = norm(halogen.iloc[i] - halogen_donor.iloc[j])
                length_bond = norm(halogen.iloc[i] - halogen_acceptor.iloc[k])
                angle = cosine(halogen.iloc[i]-halogen_donor.iloc[j], halogen.iloc[i] - halogen_acceptor.iloc[k])

                if length_bond_halogen < 2.5 and angle <-0.98 and length_bond < 3.8 and length_bond > 2.5:
                    halogen_name = ligand[ligand['X'] == halogen.iloc[i].iloc[0]].Desc.values[0]
                    halogen_donor_name = ligand[ligand['X'] == halogen_donor.iloc[j].iloc[0]].Desc.values[0]
                    halogen_acceptor_res = protein[protein['X'] == halogen_acceptor.iloc[k].iloc[0]].Res.values[0]
                    halogen_acceptor_num = protein[protein['X'] == halogen_acceptor.iloc[k].iloc[0]].Num.values[0]
                    halogen_acceptor_chain = protein[protein['X'] == halogen_acceptor.iloc[k].iloc[0]].Chain.values[0]

                    Residue.append(halogen_acceptor_res)
                    Number.append(halogen_acceptor_num)
                    Chain.append(halogen_acceptor_chain)

                    Lig_res.append(ligand[ligand['X'] == halogen.iloc[i].iloc[0]].Res.values[0])
                    Lig_num.append(ligand[ligand['X'] == halogen.iloc[i].iloc[0]].Desc.values[0])
                    Lig_chain.append(ligand[ligand['X'] == halogen.iloc[i].iloc[0]].Chain.values[0])

                    Type.append('Halogen bond')

                    distance.append(length_bond)

                    # print(f'Protein: {halogen_acceptor_res} {halogen_acceptor_num}{halogen_acceptor_chain}, Ligand: {halogen_donor_name}-{halogen_name}')
    columns = ['Residue', 'Number', 'Chain','Type']
    df_result = pd.DataFrame(columns=columns)
    df_result['Residue'] = Residue
    df_result['Number'] = Number
    df_result['Chain'] = Chain
    df_result['Type'] = Type
    df_result['Ligand'] = Lig_res
    df_result['Atom'] = Lig_num
    df_result['Distance'] = distance
    df_result = df_result.drop_duplicates()
    df_result = df_result.reset_index(drop=True, inplace=False)

    columns_ = ['Ligand', 'Atom', 'Chain','Type']
    df_pharmacophore = pd.DataFrame(columns=columns_)
    df_pharmacophore['Ligand'] = Lig_res
    df_pharmacophore['Atom'] = Lig_num
    df_pharmacophore['Chain'] = Lig_chain
    df_pharmacophore['Type'] = Type
    df_pharmacophore = df_pharmacophore.drop_duplicates()
    df_pharmacophore = df_pharmacophore.reset_index(drop=True, inplace=False)

    return(df_result, df_pharmacophore)

# print(halogen_bond_interaction(df3, df1))
