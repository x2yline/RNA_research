import os

os.chdir(r'C:\Users\Administrator\Desktop\桌面\mature.fa')
file_path_mir = 'mature.fa'

def reverse_align(sequence):
    sequence = sequence.upper().replace('U','T')
    sub_matrix = [('T','A'),('G','C'),('A','T'),('C','G')]
    for sub_item in sub_matrix:
        sequence = sequence.replace(sub_item[0],sub_item[1])
    sequence = sequence[-1::-1]
    return(sequence)

def extract_hsa_mir(file_path):
    all_miRNAs = {}
    with open('mature.fa', 'r') as f:
        find = 0
        for line in f:
            if line.startswith('>hsa'):
                miRNA_name = line.split()[0][1:]
                find = 1
            elif find == 1:
                all_miRNAs[miRNA_name] = line.strip()
            else:
                find = 0
    return(all_miRNAs)

all_miRNAs = extract_hsa_mir(file_path_mir)

def find_motif_gene(all_miRNAs,motifs=None):
    '''AGUGUU at the 3′ end of human miR-29b
    This was reflected in a patent listing UGUGUU, ACUGUU, AGAGUU, AGUCUU, 
    AGUGAU, AGUGUA, AGNGUN (Hwang et al. 2007b); notably ASUS (S = C or G) 
    appeared in six of the eight (underlined) motifs.'''
    if not motifs:
        motifs = ['AGUGUU','UGUGUU','ACUGUU','AGAGUU','AGUCUU','AGUGAU','AGUGUA']
    find_gene = {}
    for keys in all_miRNAs.keys():
        for domain in motifs:
            if all_miRNAs[keys].find(domain) != -1:
                find_gene[keys] = all_miRNAs[keys]
            elif all_miRNAs[keys].find('AGNGU')!= -1 and all_miRNAs[keys].find('AGNGU')+5 !=len(all_miRNAs[keys]):
                find_gene[keys] = all_miRNAs[keys]
    return(find_gene)
    
nuclear_miRNA = find_motif_gene(all_miRNAs)

def extract_promoter(file_path='Hs_EPDnew.dat'):
    promoter_dict = {}
    with open('Hs_EPDnew.dat') as f:
        get = 0
        for line in f:
            if line.startswith("GN"):
                gene_name = line.split('=')[-1].split(';')[0]
                if gene_name:# in gene_list:
                    get = 1
                else:
                    get = 0
            elif get == 1 and line.startswith("SE"):
                promoter_dict[gene_name] = line.split()[-1]
                get = 0
    return(promoter_dict)

promoter_dict = extract_promoter()

def alignment_miRNAs_sequences(miRNAs, sequences,seed_length=8):
    target_pairs = []
    for key1 in miRNAs.keys():
        for key2 in sequences.keys():
            if reverse_align(miRNAs[key1][1:seed_length+1]) in sequences[key2].upper().replace('U','T'):
                print('Get %s and %s\n'%(key1,key2))
                target_pairs.append((key1,key2))
    return(target_pairs)

target_pairs = alignment_miRNAs_sequences(all_miRNAs, promoter_dict,seed_length=18)
