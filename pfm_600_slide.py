'''
#----------------选出长度大于11的矩阵-------------

def parse_matrix_line(line):
    numbers_str = line.split("[")[1].split("]")[0]
    numbers = list(map(int, numbers_str.split()))
    return numbers

def read_and_filter_pwm(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    filtered_matrices = []
    current_matrix = []
    for line in lines:
        if line.startswith('>'):  
            if current_matrix:
                column_count = len(parse_matrix_line(current_matrix[1]))
                print(f"Current matrix column count: {column_count}")  
                if column_count >= 11:
                    filtered_matrices.extend(current_matrix)
                    filtered_matrices.append('')  
            current_matrix = [line.strip()]  
        elif line.strip():  
            current_matrix.append(line.strip())
    if current_matrix:
        column_count = len(parse_matrix_line(current_matrix[1]))
        if column_count >= 12:
            filtered_matrices.extend(current_matrix)

    print(f"Number of filtered matrices: {len(filtered_matrices) // 2}")  
    with open(output_file, 'w') as f:
        for line in filtered_matrices:
            f.write(f"{line}\n")


input_file_path = 'D:/tss_atac/mofit_PWM_all.txt' 
output_file_path = 'D:/tss_atac/mofit_PWM_select_1.txt'  

read_and_filter_pwm(input_file_path, output_file_path)

print("Filtered PWMs have been saved to:", output_file_path)
'''
'''
#--------选择600个矩阵----------

def select_first_n_matrices(input_file, output_file, n=600):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    selected_matrices = []
    current_matrix = []
    for line in lines:
        if line.startswith('>'):  
            if current_matrix:  
                selected_matrices.append(current_matrix)
                if len(selected_matrices) == n:  
                    break
                current_matrix = [line.strip()]  
            else:  
                current_matrix = [line.strip()]
        else:
            if line.strip(): 
                current_matrix.append(line.strip())

    if current_matrix and len(selected_matrices) < n:
        selected_matrices.append(current_matrix)
    with open(output_file, 'w') as f:
        for matrix in selected_matrices:
            for line in matrix:
                f.write(f"{line}\n")
            f.write("\n")  


# 文件的路径
input_file_path = "D:/tss_atac/mofit_PWM_select_1.txt"
output_file_path = 'D:/tss_atac/mofit_PWM_600.txt'

# 选择前600个矩阵并保存到新文件
select_first_n_matrices(input_file_path, output_file_path, n=600)

print("Selected matrices have been saved to:", output_file_path)
'''
'''
#-----------使矩阵的长度为12------------

def process_matrices(input_file, output_file, target_columns=12):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    processed_matrices = []
    current_matrix = []
    for line in lines:
        if line.startswith('>'):  
            if current_matrix:  
                processed_matrices.extend(process_matrix(current_matrix, target_columns))
                processed_matrices.append('')  
            current_matrix = [line.strip()]  
        elif line.strip():  
            current_matrix.append(line.strip())
    if current_matrix:
        processed_matrices.extend(process_matrix(current_matrix, target_columns))
    with open(output_file, 'w') as f:
        for line in processed_matrices:
            f.write(f"{line}\n")


def process_matrix(matrix, target_columns):
    #截取或补零
    processed_matrix = [matrix[0]]  
    for line in matrix[1:]:
        
        numbers = line.split('[')[-1].split(']')[0].split()
        numbers = [int(num) for num in numbers]
        adjusted_numbers = numbers[:target_columns] + [0] * (target_columns - len(numbers))
        new_line = f"  [ {' '.join(map(str, adjusted_numbers))} ]"
        processed_matrix.append(new_line)
    return processed_matrix



input_file_path = "D:/tss_atac/mofit_PWM_600.txt"
output_file_path = "D:/tss_atac/mofit_PWM_600_12.txt"

# 执行函数
process_matrices(input_file_path, output_file_path, target_columns=12)

print("Processed matrices have been saved to:", output_file_path)
'''
'''
#----------将PWM转换为PFM--------------

def normalize_matrix(matrix):
    normalized_matrix = []
    # 计算每列的总和
    column_sums = [sum(column) for column in zip(*matrix)]
    for row in matrix:
        normalized_row = [value / column_sum if column_sum else 0 for value, column_sum in zip(row, column_sums)]
        normalized_matrix.append(normalized_row)
    return normalized_matrix


def process_matrices(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    current_matrix = []
    with open(output_file, 'w') as f:
        for line in lines:
            if line.startswith('>'):  
                if current_matrix:
                    normalized_matrix = normalize_matrix(current_matrix)
                    for normalized_row in normalized_matrix:
                        f.write("  [ " + ' '.join(f"{value:.2f}" for value in normalized_row) + " ]\n")
                    f.write("\n")  
                f.write(line)  
                current_matrix = []  
            elif line.strip():  
                
                numbers = [float(num) for num in line.split('[')[-1].split(']')[0].split()]
                current_matrix.append(numbers)

        
        if current_matrix:
            normalized_matrix = normalize_matrix(current_matrix)
            for normalized_row in normalized_matrix:
                f.write("  [ " + ' '.join(f"{value:.2f}" for value in normalized_row) + " ]\n")



# 输入和输出文件路径
input_file_path = "D:/tss_atac/mofit_PWM_600_12.txt"
output_file_path = "D:/tss_atac/mofit_PFM_600_12.txt"

process_matrices(input_file_path, output_file_path)

'''

import pandas as pd
import numpy as np

#----------------计算点乘-------------------


def scan_sequences_for_motifs(promoter_one_hot_sequences, pfms):
    num_sequences = promoter_one_hot_sequences.shape[0]
    sequence_length = promoter_one_hot_sequences.shape[2]
    num_motifs = pfms.shape[0]
    motif_length = pfms.shape[2]
    num_positions = sequence_length - motif_length + 1
    promoter_match_scores = np.zeros((num_sequences, num_motifs, num_positions))

    # 对每个序列进行扫描
    for i in range(num_sequences):
        for j in range(num_motifs):
            for position in range(num_positions):
                promoter_match_scores[i, j, position] = np.sum(
                    promoter_one_hot_sequences[i, :, position:position + motif_length] * pfms[j]
                )

    return promoter_match_scores


#-----------------------------------------------------------------------------------------
#promoter

file_path = "D:/tss_atac/model_process_scan/promoter_100_one_hot_padded.csv"
# 读取文件
promoter_one_hot = pd.read_csv(file_path, header=None)
promoter_one_hot_array = promoter_one_hot.to_numpy()

num_sequences = 10000  
sequence_length = 100  
promoter_one_hot_sequences = promoter_one_hot_array.reshape(num_sequences, 4, sequence_length)

#-----------------------------------------------------------------------------------------
num_motifs = 600  
motif_length = 12  
pfms = np.zeros((num_motifs, 4, motif_length)) 
motif_index = 0  
with open('D:/tss_atac/model_process_scan/mofit_PFM_600_12.txt', 'r') as file:
    lines = file.readlines()  
    for line in lines:
        line = line.strip()  
        if line.startswith('>') or not line: 
            continue  
        cleaned_line = line.replace('[', '').replace(']', '')
        try:
            frequencies = np.array(cleaned_line.split(), dtype=float)  
            pfms[motif_index // 4, motif_index % 4, :] = frequencies  
            motif_index += 1  
        except ValueError as e:
            print(f"Error processing line: {line}\n{e}")


#-----------------------------------------------------------------------------------------
# 执行扫描
promoter_match_scores = scan_sequences_for_motifs(promoter_one_hot_sequences, pfms)
# 保存
np.save("D:/tss_atac/model_process_scan/promoter_match_scores.npy", promoter_match_scores)
np.savetxt('D:/tss_atac/model_process_scan/promoter_match_scores.csv', promoter_match_scores.reshape(-1, promoter_match_scores.shape[-1]), delimiter=',')
#---------------------------------------------------------------------------------------------


'''
#unpromoter

def scan_sequences_for_motifs(unpromoter_one_hot_sequences, pfms):
    num_sequences = unpromoter_one_hot_sequences.shape[0]
    sequence_length = unpromoter_one_hot_sequences.shape[2]
    num_motifs = pfms.shape[0]
    motif_length = pfms.shape[2]
    num_positions = sequence_length - motif_length + 1
    unpromoter_match_scores = np.zeros((num_sequences, num_motifs, num_positions))

    # 对每个序列进行扫描
    for i in range(num_sequences):
        for j in range(num_motifs):
            for position in range(num_positions):
                unpromoter_match_scores[i, j, position] = np.sum(
                    unpromoter_one_hot_sequences[i, :, position:position + motif_length] * pfms[j]
                )

    return unpromoter_match_scores

#---------------------------------------------------------------------------------------------
#unpromoter

# 文件路径
file_path = "D:/tss_atac/model_process_scan/unpromoter_100_one_hot_padded.csv"

#读取文件
unpromoter_one_hot = pd.read_csv(file_path, header=None)
unpromoter_one_hot_array = unpromoter_one_hot.to_numpy()
num_sequences = 10000  
sequence_length = 100  
unpromoter_one_hot_sequences = unpromoter_one_hot_array.reshape(num_sequences, 4, sequence_length)
#-----------------------------------------------------------------------------------------
num_motifs = 600  
motif_length = 12  
pfms = np.zeros((num_motifs, 4, motif_length))  
motif_index = 0  
with open('D:/tss_atac/model_process_scan/mofit_PFM_600_12.txt', 'r') as file:
    lines = file.readlines()  
    for line in lines:
        line = line.strip()  
        if line.startswith('>') or not line:  
            continue  
        cleaned_line = line.replace('[', '').replace(']', '')
        try:
            frequencies = np.array(cleaned_line.split(), dtype=float)  
            pfms[motif_index // 4, motif_index % 4, :] = frequencies  
            motif_index += 1  
        except ValueError as e:
            print(f"Error processing line: {line}\n{e}")


#-----------------------------------------------------------------------------------------
# 执行扫描
unpromoter_match_scores = scan_sequences_for_motifs(unpromoter_one_hot_sequences, pfms)
np.save("D:/tss_atac/model_process_scan/unpromoter_match_scores.npy", unpromoter_match_scores)

np.savetxt('D:/tss_atac/model_process_scan/unpromoter_match_scores.csv', unpromoter_match_scores.reshape(-1, unpromoter_match_scores.shape[-1]), delimiter=',')

'''


