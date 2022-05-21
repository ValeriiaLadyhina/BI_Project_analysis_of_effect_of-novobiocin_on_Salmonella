#Make list with genes names
import pandas as pd

gene_names_list = []
with open('data/sorted_p_adj_genes.txt', 'r') as file:
    for line in file:
        gene_names_list.append(line.strip().split('|')[0])

#Make dict with gene names: locus_tag

add_info = []
with open('data/genefile.txt', 'r') as file:
    for line in file:
        add_info.append(line.strip())
key_list = list(map(lambda x: x.split(';')[0].split('=')[-1], add_info))
value_list = []
for line in add_info:
    if "old_locus_tag" in line:
        value_list.append(line.split('old_locus_tag=')[-1])
    else:
        value_list.append(line.split('locus_tag=')[-1])
# old lucus tag or locus tag if old locus tag is None
value_list_2 = []
for line in value_list:
    if len(line.split(";")) > 1:
        value_list_2.append(line.split(";")[0])
    else:
        value_list_2.append(line)
id_dict = {key_list[i] : value_list_2[i] for i in range(len(key_list))}

#Make function for Module-Gene transforming

def module_transform(module_file:str):
    module = []
    with open(module_file, "r") as file:
        for line in file:
            module.append(line.strip())
    module = module[1:]
    module = list(map(lambda x: x.split('\t')[-1], module))
    module = list(map(lambda x: x.split('|')[0], module))
    new_module = list(map(lambda x: id_dict[x], module))
    with open("New_"+module_file, "w") as file:
        for line in new_module:
            file.write(line+"\n")

def rna_gene(file):
    print(file)
    new_file=pd.read_table(file)
    rna = new_file['x'].str.startswith('rna')
    for i in range(len(rna)):
        if rna[i]==True:
            new_file['x'][i]='gene'+new_file['x'][i][3::]
    new_file.to_csv('rna_'+ file,index=None, sep='\t')

rna_gene("yellow_result.txt")
rna_gene("black_result.txt")
rna_gene("brown_result.txt")
rna_gene("blue_result.txt")
rna_gene("turquoise_result.txt")
rna_gene("green_result.txt")
rna_gene("red_result.txt")


module_transform("rna_yellow_result.txt")
module_transform("rna_black_result.txt")
module_transform("rna_brown_result.txt")
module_transform("rna_blue_result.txt")
module_transform("rna_turquoise_result.txt")
module_transform("rna_green_result.txt")
module_transform("rna_red_result.txt")
