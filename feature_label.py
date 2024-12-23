import numpy
import numpy as np
#import torch
from pyjaspar import jaspardb
#import torch.nn as nn
from scipy.signal import convolve


one_hot_encoding = {
    'T': [0, 0, 0, 1],
    'G': [0, 0, 1, 0],
    'C': [0, 1, 0, 0],
    'A': [1, 0, 0, 0]
}

def read_fasta(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            # 跳过以'>'开头的行
            if not line.startswith('>'):
                sequences.append(line)
    return sequences

def sequence_to_one_hot(sequences):
    one_hot_list = []
    for seq in sequences:
        encoded_seq = [one_hot_encoding[base] for base in seq]
        one_hot_list.append(encoded_seq)
    return np.array(one_hot_list, dtype=np.float32)

promoter_file = './data/4_promoter.fasta'
no_promoter_file = './data/4_no-promoter.fasta'

promoter_sequences = read_fasta(promoter_file)
promoter_one_hot = sequence_to_one_hot(promoter_sequences)
promoter_label = np.array([1]*len(promoter_sequences))

no_promoter_sequences = read_fasta(no_promoter_file)
no_promoter_one_hot = sequence_to_one_hot(no_promoter_sequences)
no_promoter_label = np.array([0]*len(no_promoter_sequences))

data=np.concatenate((promoter_one_hot, no_promoter_one_hot), axis=0)
labels = np.concatenate((promoter_label, no_promoter_label), axis=0)
np.save('./data/labels.npy',labels)


jdb_obj = jaspardb()
motifs = jdb_obj.fetch_motifs(collection = 'CORE',tax_group = ['vertebrates'], species= 9606)
filtered_motifs = [motif for motif in motifs if motif.length == 10]
#print(len(filtered_motifs))


conv_results = []

for sample_id in range(len(data)):
    print(sample_id)
    sample = data[sample_id]
    sample_conv_result = []
    for motif in filtered_motifs:
        conv_result = np.array([
            convolve(sample[:, i], motif.pwm[i], mode='valid')
            for i in range(4)  # 针对四个碱基 A, C, G, T 单独做卷积
        ])
        sample_conv_result.append(np.sum(conv_result, axis=0).tolist())
    conv_results.append(numpy.array(sample_conv_result))

conv_results = np.array(conv_results)
np.save('./data/conv_results.npy', conv_results)


print('over')
