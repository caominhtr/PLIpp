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

def hydrogen_bond_water(water):
    acceptor = water[water['Desc'] == "O"]
    hydrogen = water[water['Name'] == "H"]
    donor = acceptor

    return donor, hydrogen, acceptor

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

def water_bridge_interaction(protein, ligand, water):
    data_protein = hydrogen_bond_protein(protein)
    data_ligand = hydrogen_bond_ligand(ligand)
    data_water = hydrogen_bond_water(water)


    protein_donor = data_protein[0]
    protein_hydrogen = data_protein[1]
    protein_acceptor = data_protein[2]

    ligand_donor = data_ligand[0]
    ligand_hydrogen = data_ligand[1]
    ligand_acceptor = data_ligand[2]

    water_donor = data_water[0]
    water_hydrogen = data_water[1]
    water_acceptor = data_water[2]

    D1 = []     #protein donor
    H1 = []     #protein H   
    A1 = np.array(protein_acceptor[['X','Y','Z']])
    D2 = []     #water donor
    H2 = []     #water H
    A2 = np.array(water_acceptor[['X','Y','Z']])
    D3 = []     #ligand donor
    H3 = []     #ligand H
    A3 = np.array(ligand_acceptor[['X','Y','Z']])

    for order1 in protein_donor.Order:
        num1 = protein_donor[protein_donor['Order'] == order1]['Num'].values[0]
        for order2 in protein_hydrogen[protein_hydrogen['Num']== num1].Order:
            D_1 = protein_donor[(protein_donor['Order']==order1) & (protein_donor['Num']==num1)][['X','Y','Z']].values[0]
            H_1 = protein_hydrogen[(protein_hydrogen['Num']==num1)&(protein_hydrogen['Order']==order2)][['X','Y','Z']].values[0]
            if (norm(D_1-H_1)<1.1):
                D1.append(D_1)
                H1.append(H_1)

    for order3 in water_donor.Order:
        num2 = water_donor[water_donor['Order'] == order3]['Num'].values[0]
        for order4 in water_hydrogen[water_hydrogen['Num']== num2].Order:
            D_2 = water_donor[(water_donor['Order']==order3) & (water_donor['Num']==num2)][['X','Y','Z']].values[0]
            H_2 = water_hydrogen[(water_hydrogen['Num']==num2)&(water_hydrogen['Order']==order4)][['X','Y','Z']].values[0]
            if (norm(D_2-H_2)<1.1):
                D2.append(D_2)
                H2.append(H_2)

    for order4 in ligand_donor.Order:
        num3 = ligand_donor[ligand_donor['Order'] == order4]['Num'].values[0]
        for order5 in ligand_hydrogen[ligand_hydrogen['Num']== num3].Order:
            D_3 = ligand_donor[(ligand_donor['Order']==order4) & (ligand_donor['Num']==num3)][['X','Y','Z']].values[0]
            H_3 = ligand_hydrogen[(ligand_hydrogen['Num']==num3)&(ligand_hydrogen['Order']==order5)][['X','Y','Z']].values[0]
            if (norm(D_3-H_3)<1.1):
                D3.append(D_3)
                H3.append(H_3)
    
    Residue = []
    Number = []
    Chain = []
    Type = []

    Lig_res = []
    Lig_desc = []
    Lig_chain = []
    Lig_type = []

    distance = []

    # Case 1: Protein (D_1) - Protein(H_1) - Water (D_2) - Water(H_2) - Ligand (A_3)
    for i in range(len(D1)):
        for j in range(len(D2)):
            for k in range(len(A3)):
                if np.logical_and(np.logical_and((norm(D1[i]-D2[j])>2.5),(norm(D1[i]-D2[j])<3.8)), np.logical_and((norm(D2[j]-A3[k])>2.5),(norm(D2[j]-A3[k])<3.8))):
                    if np.logical_and(cosine(H1[i]-D1[i], H1[i]-D2[j]) < -0.75, cosine(H2[j]-D2[j], H2[j]-A3[k])<-0.65):
                        Prores = protein[protein['X'] == D1[i][0]].Res.values[0]
                        Pronum = protein[protein['X'] == D1[i][0]].Num.values[0]
                        Prochain = protein[protein['X'] == D1[i][0]].Chain.values[0]
                        Waternum = water[water['X'] == D2[j][0]].Num.values[0]
                        liganddesc = ligand[ligand['X'] == A3[k][0]].Desc.values[0]
                        ligandchain = ligand[ligand['X'] == A3[k][0]].Chain.values[0]

                        Residue.append(Prores)
                        Number.append(Pronum)
                        Chain.append(Prochain)
                        Type.append('Water bridge donor')

                        Lig_res.append(ligand[ligand['X'] == A3[k][0]].Res.values[0])
                        Lig_desc.append(liganddesc)
                        Lig_chain.append(ligandchain)
                        Lig_type.append('Water bridge acceptor')

                        distance.append([norm(D1[i]-D2[j]), norm(D2[j]-A3[k])])

                        
                        # print(f'Protein (Donor): {Prores} {Pronum}{Prochain}, Ligand (Acceptor): {liganddesc} chain {ligandchain}, Water (Donor): HOH{Waternum}')

    # Case 2: Protein (D_1) - Protein(H_1) - Water (A_2) - Water(H_3) - Ligand (D_3)
    for i in range(len(D1)):
        for j in range(len(D3)):
            for k in range(len(A2)):
                if np.logical_and(np.logical_and((norm(D1[i]-A2[k])>2.5),(norm(D1[i]-A2[k])<3.8)), np.logical_and((norm(D3[j]-A2[k])>2.5),(norm(D3[j]-A2[k])<3.8))):
                    if np.logical_and(cosine(H1[i]-D1[i], H1[i]-A2[k]) < -0.75, cosine(H3[j]-A2[k], H3[j]-D3[j])<-0.65):
                        Prores = protein[protein['X'] == D1[i][0]].Res.values[0]
                        Pronum = protein[protein['X'] == D1[i][0]].Num.values[0]
                        Prochain = protein[protein['X'] == D1[i][0]].Chain.values[0]
                        Waternum = water[water['X'] == A2[k][0]].Num.values[0]
                        liganddesc = ligand[ligand['X'] == D3[j][0]].Desc.values[0]
                        ligandchain = ligand[ligand['X'] == D3[j][0]].Chain.values[0]


                        Residue.append(Prores)
                        Number.append(Pronum)
                        Chain.append(Prochain)
                        Type.append('Water bridge donor')

                        Lig_res.append(ligand[ligand['X'] == D3[j][0]].Res.values[0])
                        Lig_desc.append(liganddesc)
                        Lig_chain.append(ligandchain)
                        Lig_type.append('Water bridge donor')

                        distance.append([norm(D1[i]-A2[k]), norm(D3[j]-A2[k])])

                        # print(f'Protein (Donor): {Prores} {Pronum}{Prochain}, Ligand (Donor): {liganddesc} chain {ligandchain}, Water (Acceptor): HOH{Waternum}')


    # Case 3: Protein (A_1) - Water(H_2) - Water (D_2) - Ligand (H_3) - Ligand (D_3)
    for i in range(len(A1)):
        for j in range(len(D3)):
            for k in range(len(D2)):
                if np.logical_and(np.logical_and((norm(A1[i]-D2[k])>2.5),(norm(A1[i]-D2[k])<3.8)), np.logical_and((norm(D3[j]-D2[k])>2.5),(norm(D3[j]-D2[k])<3.8))):
                    if np.logical_and(cosine(H2[k]-A1[i], H2[k]-D2[k]) < -0.75, cosine(H3[j]-D2[k], H3[j]-D3[j])<-0.65):
                        Prores = protein[protein['X'] == A1[i][0]].Res.values[0]
                        Pronum = protein[protein['X'] == A1[i][0]].Num.values[0]
                        Prochain = protein[protein['X'] == A1[i][0]].Chain.values[0]
                        Waternum = water[water['X'] == D2[k][0]].Num.values[0]
                        liganddesc = ligand[ligand['X'] == D3[j][0]].Desc.values[0]
                        ligandchain = ligand[ligand['X'] == D3[j][0]].Chain.values[0]


                        Residue.append(Prores)
                        Number.append(Pronum)
                        Chain.append(Prochain)
                        Type.append('Water bridge acceptor')

                        Lig_res.append(ligand[ligand['X'] == D3[j][0]].Res.values[0])
                        Lig_desc.append(liganddesc)
                        Lig_chain.append(ligandchain)
                        Lig_type.append('Water bridge donor')

                        distance.append([norm(A1[i]-D2[k]), norm(D3[j]-D2[k])])
                        # print(f'Protein (Acceptor): {Prores} {Pronum}{Prochain}, Ligand (Donor): {liganddesc} chain {ligandchain}, Water (Donor): HOH{Waternum}')

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

    return(df_result, df_pharmacophore)

# print(water_bridge_interaction(df3, df1, df5)[0])
# print(water_bridge_interaction(df3, df1, df5)[1])
