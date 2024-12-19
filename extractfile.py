import sys
import os
import shutil
import pandas as pd
import numpy as np

if len(sys.argv) < 2:
    print("Please provide a file path as an argument.")
    sys.exit(1)

file_path = sys.argv[1]

base_name = os.path.splitext(os.path.basename(file_path))[0]

output_directory = base_name
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

try:
    shutil.copy(file_path, os.path.join(output_directory, os.path.basename(file_path)))

except FileNotFoundError:
    print(f"Error: The input file {file_path} does not exist.")
    sys.exit(1)
except Exception as e:
    print(f"An unexpected error occurred: {e}")
    sys.exit(1)

os.chdir(output_directory)


def make_protein_file(input_file_path):
    try:
        output_protein_file_path = "protein.txt"  
        with open(output_protein_file_path, 'w') as output_file:
            with open(input_file_path, 'r') as input_file:
                for line in input_file:
                    if line.startswith("ATOM"):
                        output_file.write(line) 

        return output_protein_file_path  

    except FileNotFoundError:
        print(f"Error: The input file {input_file_path} does not exist.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None
    
def make_ligand_file(input_file_path):
    try:
        output_ligand_file_path = "ligand.txt"  
        with open(output_ligand_file_path, 'w') as output_file:
            with open(input_file_path, 'r') as input_file:
                for line in input_file:
            
                    if line.startswith("HETATM") and "HOH" not in line:
                        output_file.write(line)  
        return output_ligand_file_path  

    except FileNotFoundError:
        print(f"Error: The input file {input_file_path} does not exist.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None
    
def make_water_file(input_file_path):
    try:
        output_water_file_path = "water.txt"  
        with open(output_water_file_path, 'w') as output_file:
            with open(input_file_path, 'r') as input_file:
                for line in input_file:
            
                    if line.startswith("HETATM") and "HOH" in line:
                        output_file.write(line)  
        return output_water_file_path  

    except FileNotFoundError:
        print(f"Error: The input file {input_file_path} does not exist.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None
    
def make_conect_file(input_file_path):
    try:
        output_conect_file_path = "conect.txt"  
        with open(output_conect_file_path, 'w') as output_file:
            with open(input_file_path, 'r') as input_file:
                for line in input_file:
            
                    if line.startswith("CONECT"):
                        output_file.write(line)  
        return output_conect_file_path  

    except FileNotFoundError:
        print(f"Error: The input file {input_file_path} does not exist.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None
def read_csv_with_null(file_path, **kwargs):
    try:
        df = pd.read_csv(file_path, **kwargs)
    except pd.errors.EmptyDataError:
        df = None
    return df

def process_file(file_path):
    file_path_ = f"{base_name}_filtered.txt"
    with open(file_path, 'r') as infile:
        lines = infile.readlines()

    filtered_lines = [line for line in lines if 'ANISOU' not in line]

    with open(file_path_, 'w') as outfile:
        outfile.writelines(filtered_lines)


    output_protein_file_path = make_protein_file(file_path_)
    output_ligand_file_path = make_ligand_file(file_path_)
    output_water_file_path = make_water_file(file_path_)
    output_conect_file_path = make_conect_file(file_path_)

    if None in [output_protein_file_path, output_ligand_file_path, output_water_file_path, output_conect_file_path]:
        print("Error: One or more required files were not created.")
        return None 
    
    try:
        df = read_csv_with_null('ligand.txt', sep=r'\s+', header=None)
        df2 = read_csv_with_null('protein.txt', sep=r'\s+', header=None)
        df4 = read_csv_with_null('water.txt', sep=r'\s+', header=None)

    except Exception as e:
        print(f"Error reading one of the required files: {e}")
        return None

    list_check = []
    for i in range(df.shape[0]):
        if any(df.iloc[i].isna()):
            list_check.append(df.iloc[i][1])
    for j in range(df2.shape[0]):
        if any(df2.iloc[j].isna()):
            list_check.append(df2.iloc[j][1])
    if df4 is not None:
        for k in range(df4.shape[0]):
            if any(df4.iloc[k].isna()):
                list_check.append(df4.iloc[k][1])
    if len(list_check) > 0:
        print(f"Check the pdb file at the atom {list_check}")
    if df.shape[1] != 12:
        split_col = df.iloc[:, 4].str.extract(r'([A-Z])(\d+)')
        df.insert(5, 'Num', split_col[1])
        df.iloc[:, 4] = split_col[0]
        col_to_drop_ = [0, 8, 9]
        df1 = df.drop(columns= col_to_drop_, axis = 1)
        df1.columns = ['Order','Desc', 'Res', 'Chain', 'Num', 'X', 'Y', 'Z','Name']
    else:
        col_to_dropp = [0, 9, 10]
        df1 = df.drop(columns=col_to_dropp, axis=1)
        df1.columns = ['Order', 'Desc', 'Res', 'Chain', 'Num', 'X', 'Y', 'Z', 'Name']

    if df2.shape[1] != 12:
        print(f"Check the protein file")
    if df4 is not None:
        if df4.shape[1] != 12:
            split_col = df4.iloc[:, 4].str.extract(r'([A-Z])(\d+)')
            df4.insert(5, 'Num', split_col[1])
            df4.iloc[:, 4] = split_col[0]
            col_to_drop_ = [0, 8, 9]
            df5 = df4.drop(columns= col_to_drop_, axis = 1)
            df5.columns = ['Order','Desc', 'Res', 'Chain', 'Num', 'X', 'Y', 'Z','Name']
        else:
            col_to_dropp = [0, 9, 10]
            df5 = df4.drop(columns=col_to_dropp, axis=1)
            df5.columns = ['Order', 'Desc', 'Res', 'Chain', 'Num', 'X', 'Y', 'Z', 'Name']
            
    col_to_drop = [0, 9, 10]
    df3 = df2.drop(columns=col_to_drop, axis=1)
    df3.columns = ['Order', 'Desc', 'Res', 'Chain', 'Num', 'X', 'Y', 'Z', 'Name']
    
    if df4 is None:
        df5 = pd.DataFrame(columns= ['Order', 'Desc', 'Res', 'Chain', 'Num', 'X', 'Y', 'Z', 'Name'])
    

    if df1.Name.isna().any():
        return None
    elif df3.Name.isna().any():
        return None

    data_list = []
    list_new_indices = []
    try:
        with open('conect.txt', 'r') as file:
            for line in file:
                parts = line.strip().split()
                if len(parts) > 1:
                    numbers_list = [int(part) for part in parts[1:]]
                    data_list.append(numbers_list)
    except Exception as e:
        print(f"Error reading conect file: {e}")
        return None

    list_new = []
    for index in data_list:
        if not any(i not in df1['Order'].values for i in index):
            custom_order = {val: idx for idx, val in enumerate(index)}
            list_new.append(np.array(df1[df1['Order'].isin(custom_order)].sort_values(by='Order', key=lambda x: x.map(custom_order))['Name']))
            list_new_indices.append(data_list.index(index))
    
 

    return df1, df3, df5, data_list, list_new_indices, list_new



df1, df3, df5, data_list, list_new_indices, list_new = process_file(file_path)           



