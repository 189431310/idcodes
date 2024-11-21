# -*- coding: utf-8 -*-

from numpy.random import seed
seed(8) 

import tensorflow
tensorflow.random.set_seed(7)


import numpy as np 
import pandas as pd 
import tensorflow as tf
import os

from tensorflow.keras import backend as K
from tensorflow.keras.models import Model ,load_model
from tensorflow.keras.layers import Flatten, Dense, Dropout
from tensorflow.keras.applications.inception_resnet_v2 import InceptionResNetV2, preprocess_input
from keras.applications.vgg16 import preprocess_input
from keras.applications.vgg16 import decode_predictions
from keras.applications.vgg16 import VGG16
from tensorflow.keras.optimizers import Adam, RMSprop
from tensorflow.keras.preprocessing.image import ImageDataGenerator
from tensorflow.keras.callbacks import ModelCheckpoint
import numpy as np


from tensorflow.python.keras import models
from tensorflow.python.keras import layers
from tensorflow.keras import optimizers


import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

#---------------------------------------------------------------------------------------------------------------------
from scipy.ndimage import rotate

#data
positive_samples = pd.read_csv("/home/1015_1/project/jyli/promoter_match_scores.csv", header=None)
negative_samples = pd.read_csv("/home/1015_1/project/jyli/unpromoter_match_scores.csv", header=None)


#添加索引列
positive_samples['index'] = positive_samples.index
negative_samples['index'] = negative_samples.index


positive_indices = positive_samples['index'].values
negative_indices = negative_samples['index'].values
positive_samples.drop(columns='index', inplace=True)
negative_samples.drop(columns='index', inplace=True)

# 调整数据形状
positive_reshaped = positive_samples.values.reshape((-1, 600, 89))
negative_reshaped = negative_samples.values.reshape((-1, 600, 89))


#正样本索引
positive_indices = np.arange(1,len(positive_reshaped)+1)
#负样本索引
negative_indices = np.arange(len(positive_reshaped)+1, len(positive_reshaped) + len(negative_reshaped)+1)

# 合并索引
indices = np.concatenate((positive_indices, negative_indices), axis=0)

# 定义旋转角度范围
angle_range = (-10, 10)


# 对正样本进行随机旋转
for i in range(positive_reshaped.shape[0]):
    angle = np.random.uniform(angle_range[0], angle_range[1])  # 随机生成旋转角度
    positive_reshaped[i] = rotate(positive_reshaped[i], angle, reshape=False)

# 对负样本进行随机旋转
for i in range(negative_reshaped.shape[0]):
    angle = np.random.uniform(angle_range[0], angle_range[1])
    negative_reshaped[i] = rotate(negative_reshaped[i], angle, reshape=False)

# 打印数据形状
print("Positive data shape:", positive_reshaped.shape)
print("Negative data shape:", negative_reshaped.shape)

#------------------------------------------------------------------------

import numpy as np

def apply_random_mask(matrix, block_size=(30, 10), fill_value=0):
    H, W = matrix.shape
    h, w = block_size
    h = min(h, H)
    w = min(w, W)
    start_i = np.random.randint(0, H - h + 1)
    start_j = np.random.randint(0, W - w + 1)

    # 创建遮挡区域
    matrix[start_i:start_i + h, start_j:start_j + w] = fill_value

    return matrix

#随机遮挡到正样本和负样本
for i in range(positive_reshaped.shape[0]):
    positive_reshaped[i] = apply_random_mask(positive_reshaped[i], block_size=(30, 10), fill_value=0)

for i in range(negative_reshaped.shape[0]):
    negative_reshaped[i] = apply_random_mask(negative_reshaped[i], block_size=(30, 10), fill_value=0)

# 打印数据形状
print("Positive data random_ shape:", positive_reshaped.shape)
print("Negative data random_ shape:", negative_reshaped.shape)



#------------------------------------------------------------------------


positive_labels = np.ones((positive_reshaped.shape[0],), dtype=int)
negative_labels = np.zeros((negative_reshaped.shape[0],), dtype=int)

X_po_train, X_po_temp, y_po_train, y_po_temp, indices_po_train, indices_po_temp = train_test_split(positive_reshaped, positive_labels, positive_indices, test_size=0.4, random_state=50)
X_po_test, X_po_validation, y_po_test, y_po_validation, indices_po_test, indices_po_validation = train_test_split(X_po_temp, y_po_temp, indices_po_temp, test_size=0.5, random_state=50)

X_ne_train, X_ne_temp, y_ne_train, y_ne_temp, indices_ne_train, indices_ne_temp = train_test_split(negative_reshaped, negative_labels, negative_indices, test_size=0.4, random_state=50)
X_ne_test, X_ne_validation, y_ne_test, y_ne_validation, indices_ne_test, indices_ne_validation = train_test_split(X_ne_temp, y_ne_temp, indices_ne_temp, test_size=0.5, random_state=50)

X_train = np.concatenate((X_po_train, X_ne_train), axis=0)
y_train = np.concatenate((y_po_train, y_ne_train), axis=0)
indices_train = np.concatenate((indices_po_train, indices_ne_train), axis=0)

X_test = np.concatenate((X_po_test, X_ne_test), axis=0)
y_test = np.concatenate((y_po_test, y_ne_test), axis=0)
indices_test = np.concatenate((indices_po_test, indices_ne_test), axis=0)

X_validation = np.concatenate((X_po_validation, X_ne_validation), axis=0)
y_validation = np.concatenate((y_po_validation, y_ne_validation), axis=0)
indices_validation = np.concatenate((indices_po_validation, indices_ne_validation), axis=0)


X = np.concatenate((X_train, X_test, X_validation), axis=0)
y = np.concatenate((y_train, y_test, y_validation), axis=0)
indices = np.concatenate((indices_train, indices_test, indices_validation), axis=0)

print("X shape:", X.shape)
print("y shape:", y.shape)

test_indices_df = pd.DataFrame(indices_test, columns=['Test Indices'])
validation_indices_df = pd.DataFrame(indices_validation, columns=['Validation Indices'])

# 保存
test_indices_df.to_csv('test__indices.csv', index=False)
validation_indices_df.to_csv('validation__indices.csv', index=False)

print("?????????:", X_train.shape, y_train.shape)
print("?????????:", X_test.shape, y_test.shape)
print("?????????:", X_validation.shape, y_validation.shape)



batch_size = 128

buffer_size = 6
train_dataset = tf.data.Dataset.from_tensor_slices((X_train, y_train))
train_dataset = train_dataset.shuffle(buffer_size).batch(batch_size)


valid_dataset = tf.data.Dataset.from_tensor_slices((X_validation, y_validation)).batch(batch_size)
test_dataset = tf.data.Dataset.from_tensor_slices((X_test, y_test)).batch(batch_size)


#-------------------------------------------------------------------------------------



num_train_batches = tf.data.experimental.cardinality(train_dataset).numpy()
num_valid_batches = tf.data.experimental.cardinality(valid_dataset).numpy()



from tensorflow.keras.layers import Dense, Dropout, Conv1D, MaxPooling1D, Flatten
from tensorflow.keras import models, layers
from keras import regularizers
from tensorflow.keras.regularizers import l2, l1
from keras.layers import Conv2D

#---------------------------------------------------------------------------------------------
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Add, Conv2D, MaxPooling2D, Flatten, Dense, Dropout,BatchNormalization, Activation
from tensorflow.keras.optimizers import SGD
from tensorflow.keras.initializers import Zeros
from tensorflow.keras.initializers import RandomNormal

from tensorflow.keras.layers import AveragePooling2D
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, ReLU, LeakyReLU
from tensorflow.keras.optimizers import Adam



# 定义输入
input_tensor = Input(shape=(600, 89, 1))

# 第一个卷积层
x = Conv2D(32, (7, 7), strides=(1, 1), padding='same', kernel_regularizer=l2(0.0025))(input_tensor)
x = BatchNormalization()(x)
x = LeakyReLU()(x)

# 保存第一个卷积层的输出，作为跨层连接的输入
residual = x

# 添加第二个卷积层
#x = Conv2D(32, (5, 5), padding='same', kernel_regularizer=l2(0.00096))(x)
x = Conv2D(64, (3, 3), padding='same')(x)
x = BatchNormalization()(x)
x = LeakyReLU()(x)


# 残差连接
x = Add()([residual, x])
x = LeakyReLU()(x)

# 添加池化层
x = AveragePooling2D((2, 2))(x)

# 将二维特征图展平为一维向量
x = Flatten()(x)
x = Dropout(0.3)(x)
x = Dense(256)(x)
x = BatchNormalization()(x)
x = LeakyReLU()(x)

# 添加输出层，使用Sigmoid激活函数
output_tensor = Dense(1, activation='sigmoid')(x)

# 创建模型
model = Model(inputs=input_tensor, outputs=output_tensor)

# 设置Adam优化器，学习率为0.01
optimizer = Adam(learning_rate=0.01)

# 编译模型，指定损失函数和评估指标
model.compile(loss='binary_crossentropy',
              optimizer=optimizer,
              metrics=['accuracy'])

# 打印模型结构摘要
print(model.summary())
NUM_EPOCHS=200
STEP_SIZE_TRAIN = num_train_batches
STEP_SIZE_VALID = num_valid_batches

#model.fit
result = model.fit(
    x=train_dataset,
    steps_per_epoch=STEP_SIZE_TRAIN,
    validation_data=valid_dataset,
    validation_steps=STEP_SIZE_VALID,
    epochs=NUM_EPOCHS,
    verbose=1
    )

test_loss, test_accuracy = model.evaluate(test_dataset)
print("Test Loss:", test_loss)
print("Test Accuracy:", test_accuracy)


# 验证
validation_predictions = model.predict(valid_dataset).flatten()
# 测试
test_predictions = model.predict(test_dataset).flatten()

# 结合验证集的索引和预测结果，添加标识列
validation_1results = pd.DataFrame({
    'subseq_index': indices_validation,
    'prediction': validation_predictions,
})

# 结合测试集的索引和预测结果，并添加一个标识列
test_1results = pd.DataFrame({
    'subseq_index': indices_test,
    'prediction': test_predictions,
})

# 加载每个原始序列的子序列索引的CSV文件，并处理子序列索引
orig_to_subseq_df = pd.read_csv('a_index_num_list.csv')

# 处理分号分隔的子序列索引，扩展为多行
orig_to_subseq_df['subseq_indices'] = orig_to_subseq_df['subseq_indices'].str.split('; ')
orig_to_subseq_df = orig_to_subseq_df.explode('subseq_indices').reset_index(drop=True)
orig_to_subseq_df['subseq_indices'] = orig_to_subseq_df['subseq_indices'].astype(int)

# 将预测结果数据框与原始序列到子序列的索引数据框合并
validation_merged_df = pd.merge(orig_to_subseq_df, validation_1results, left_on='subseq_indices',
                                right_on='subseq_index')
test_merged_df = pd.merge(orig_to_subseq_df, test_1results, left_on='subseq_indices', right_on='subseq_index')

# 计算原始序列的所有子序列的预测结果总和
validation_grouped_df = validation_merged_df.groupby('orig_seq_index')['prediction'].sum()
test_grouped_df = test_merged_df.groupby('orig_seq_index')['prediction'].sum()

# 判断原始序列是否为启动子
validation_grouped_df = validation_grouped_df >= 1
test_grouped_df = test_grouped_df >= 1


ori_index = pd.read_csv('ori_index.csv')
ori_label = pd.read_csv('ori_label.csv')

true_labels_df = pd.DataFrame({
    'origseqindex': [ori_index],  
    'true_label': [ori_label]  
})

# 合并预测结果和真实标签
validation_merged_true = pd.merge(validation_grouped_df, true_labels_df, left_index=True, right_on='origseqindex')
test_merged_true = pd.merge(test_grouped_df, true_labels_df, left_index=True, right_on='origseqindex')

# 计算准确率
validation_accuracy = (validation_merged_true['prediction'] == validation_merged_true['true_label']).mean()
test_accuracy = (test_merged_true['prediction'] == test_merged_true['true_label']).mean()

print(f'Validation Accuracy: {validation_accuracy}')
print(f'Test Accuracy: {test_accuracy}')

validation_predictions_df = validation_grouped_df > 0
test_predictions_df = test_grouped_df > 0

# 合并真实标签和预测结果
validation_merged_df = pd.merge(true_labels_df, validation_predictions_df, left_on='origseqindex', right_index=True)
test_merged_df = pd.merge(true_labels_df, test_predictions_df, left_on='origseqindex', right_index=True)

validation_merged_df.rename(columns={'prediction': 'predicted_label'}, inplace=True)
test_merged_df.rename(columns={'prediction': 'predicted_label'}, inplace=True)


def calculate_metrics(df):
    TP = ((df['true_label'] == 1) & (df['predicted_label'] == 1)).sum()
    TN = ((df['true_label'] == 0) & (df['predicted_label'] == 0)).sum()
    FP = ((df['true_label'] == 0) & (df['predicted_label'] == 1)).sum()
    FN = ((df['true_label'] == 1) & (df['predicted_label'] == 0)).sum()
    SN = TP / (TP + FN) if (TP + FN) > 0 else 0
    SP = TN / (TN + FP) if (TN + FP) > 0 else 0
    MCC = (TP * TN - FP * FN) / ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5 if (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) > 0 else 0
    return SN, SP, MCC
# 计算SN、SP和MCC
validation_SN, validation_SP, validation_MCC = calculate_metrics(validation_merged_df)
test_SN, test_SP, test_MCC = calculate_metrics(test_merged_df)
print(f'Validation SN: {validation_SN}, SP: {validation_SP}, MCC: {validation_MCC}')
print(f'Test SN: {test_SN}, SP: {test_SP}, MCC: {test_MCC}')



import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
test_true_labels = y_test
validation_true_labels = y_validation

fpr_test, tpr_test, thresholds_test = roc_curve(test_true_labels, test_predictions)
roc_auc_test = auc(fpr_test, tpr_test)

fpr_validation, tpr_validation, thresholds_validation = roc_curve(validation_true_labels, validation_predictions)
roc_auc_validation = auc(fpr_validation, tpr_validation)
# 绘制ROC
import matplotlib.pyplot as plt

plt.figure()
plt.plot(fpr_test, tpr_test, color='purple', lw=2, label='Test ROC curve (area = %0.2f)' % roc_auc_test)
plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')  
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
plt.legend(loc="lower right")
plt.text(1, 0, 'AUC = %0.3f' % roc_auc_test, horizontalalignment='right', verticalalignment='bottom', fontsize=12)
plt.show()


plt.figure()
plt.plot(fpr_validation, tpr_validation, color='Indigo', lw=2, label='Validation ROC curve (area = %0.2f)' % roc_auc_validation)
plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')  
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
plt.legend(loc="lower right")
plt.text(1, 0, 'AUC = %0.3f' % roc_auc_test, horizontalalignment='right', verticalalignment='bottom', fontsize=12)
plt.show()


