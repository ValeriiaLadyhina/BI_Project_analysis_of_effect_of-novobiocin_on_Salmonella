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
reverse_id_dict = {}
for key in id_dict.keys():
    reverse_id_dict[id_dict[key]] = key

def reverse_id(gene_list:list, name:str):
    with open(name, 'w') as file:
        for gene in gene_list:
            file.write(reverse_id_dict[gene]+'\n')

def get_normal_names(path):
    tcs = []
    with open(path, 'r') as file:
        for line in file:
            tcs.append(line.strip())
    tcs = tcs[1:]
    tcs = list(map(lambda x: x.split('\t')[1], tcs))
    return tcs

pathways = ["c5_branched_dibasic_acid_metabolism.txt",
            "oxocarboxylic_acid_metabolism.txt",
            "glycolysis_gluconeogenesis.txt",
            "two_component_system.txt"]
for path in pathways:
    reverse_id(get_normal_names(path), "New_"+path)
