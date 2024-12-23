import pandas as pd
from sklearn.model_selection import train_test_split
import torch
from torch.utils.data import DataLoader, TensorDataset
import torch.nn as nn
from model import *
from torch.optim.lr_scheduler import ReduceLROnPlateau
import numpy as np
from torchvision import models, transforms
from sklearn.metrics import classification_report, roc_auc_score, roc_curve, matthews_corrcoef


if __name__=='__main__':
    data= np.load('./data/conv_results.npy')
    labels = np.load('./data/labels.npy')
    test_val_results = pd.DataFrame()

    #set_seed(22)
    # 按正负样本分开
    positive_data = data[labels == 1]
    negative_data = data[labels == 0]

    # 获取正负样本的原始索引
    positive_indices = np.where(labels == 1)[0]
    negative_indices = np.where(labels == 0)[0]

    # 每部分数据的样本数
    train_size_pos = int(0.6 * len(positive_data))  # 60% 的正样本
    val_size_pos = int(0.2 * len(positive_data))  # 20% 的正样本
    test_size_pos = len(positive_data) - train_size_pos - val_size_pos  # 剩下 20%

    train_size_neg = int(0.6 * len(negative_data))  # 60% 的负样本
    val_size_neg = int(0.2 * len(negative_data))  # 20% 的负样本
    test_size_neg = len(negative_data) - train_size_neg - val_size_neg  # 剩下 20%

    # 划分正样本
    # train_pos, temp_pos = train_test_split(positive_data, train_size=train_size_pos, random_state=42)
    # val_pos, test_pos = train_test_split(temp_pos, train_size=val_size_pos, random_state=42)

    train_pos, temp_pos, train_pos_idx, temp_pos_idx = train_test_split(
        positive_data, positive_indices, train_size=train_size_pos, random_state=42
    )
    val_pos, test_pos, val_pos_idx, test_pos_idx = train_test_split(
        temp_pos, temp_pos_idx, train_size=val_size_pos, random_state=42
    )


    # 划分负样本
    # train_neg, temp_neg = train_test_split(negative_data, train_size=train_size_neg, random_state=42)
    # val_neg, test_neg = train_test_split(temp_neg, train_size=val_size_neg, random_state=42)

    train_neg, temp_neg, train_neg_idx, temp_neg_idx = train_test_split(
        negative_data, negative_indices, train_size=train_size_neg, random_state=42
    )
    val_neg, test_neg, val_neg_idx, test_neg_idx = train_test_split(
        temp_neg, temp_neg_idx, train_size=val_size_neg, random_state=42
    )

    test_val_results['val_index'] = val_pos_idx.tolist() + val_neg_idx.tolist()
    test_val_results['test_index'] = test_pos_idx.tolist() + test_neg_idx.tolist()

    # 合并正负样本
    train_data = np.concatenate((train_pos, train_neg))
    train_labels = np.array([1] * len(train_pos) + [0] * len(train_neg))

    val_data = np.concatenate((val_pos, val_neg))
    val_labels = np.array([1] * len(val_pos) + [0] * len(val_neg))

    test_data = np.concatenate((test_pos, test_neg))
    test_labels = np.array([1] * len(test_pos) + [0] * len(test_neg))

    # 打乱训练集
    train_indices = np.random.permutation(len(train_data))
    # val_indices = np.random.permutation(len(val_data))
    # test_indices = np.random.permutation(len(test_data))

    train_data, train_labels = train_data[train_indices], train_labels[train_indices]


    # val_data, val_labels = val_data[val_indices], val_labels[val_indices]
    # test_data, test_labels = test_data[test_indices], test_labels[test_indices]

    print(f"训练集样本数: {len(train_data)}, 验证集样本数: {len(val_data)}, 测试集样本数: {len(test_data)}")

# _______________________________________________转换成Tensor___________________________________
    train_data = torch.tensor(train_data, dtype=torch.float32).unsqueeze(1)  # 添加通道维度
    train_labels = torch.tensor(train_labels, dtype=torch.long)

    val_data = torch.tensor(val_data, dtype=torch.float32).unsqueeze(1)
    val_labels = torch.tensor(val_labels, dtype=torch.long)

    test_data = torch.tensor(test_data, dtype=torch.float32).unsqueeze(1)
    test_labels = torch.tensor(test_labels, dtype=torch.long)

    # 创建DataLoader
    batch_size = 64

    train_loader = DataLoader(TensorDataset(train_data, train_labels), batch_size=batch_size, shuffle=True,num_workers=8)
    val_loader = DataLoader(TensorDataset(val_data, val_labels), batch_size=batch_size, shuffle=False,num_workers=8)
    test_loader = DataLoader(TensorDataset(test_data, test_labels), batch_size=batch_size, shuffle=False,num_workers=8)

    # 打印数据维度检查
    for data, labels in train_loader:
        print("输入维度:", data.shape)  # 应该为 (batch_size, 1, 101, 72)
        print("标签维度:", labels.shape)  # 应该为 (batch_size,)
        break

    # 定义损失函数和优化器
    model =  ResNet_18(BasicBlock, num_classes=2)

#___________________________________加载已保存的模型________________________________________________________________
    model.load_state_dict(torch.load('./history_model/best_model.pth'))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    val_acc,val_loss,val_predict = evaluate(model, val_loader,device)
    test_acc,test_loss,test_predict = evaluate(model, test_loader,device)
    print(f"验证集准确率: {val_acc:.2f}%, Loss:{val_loss:.4f}")
    print(f"测试集准确率: {test_acc:.2f}%, Loss:{test_loss:.4f}")

    all_scores = {}
    for i in range(2):
        if i == 0:
            predictions = val_predict['predicted_label']
            predict_proba =val_predict['predicted_probability']
            val_labels = val_labels

            test_val_results['val_predictions']=predictions
            test_val_results['val_labels']=val_labels
        else:
            predictions = test_predict['predicted_label']
            predict_proba = test_predict['predicted_probability']
            val_labels = test_labels

            test_val_results['test_predictions']=predictions
            test_val_results['test_labels']=val_labels

        report = classification_report(val_labels, predictions, digits=5)
        print(report)
        # sn = report["0"]["recall"]
        # sp = report["1"]["recall"]
        # acc = report['accuracy']
        # MCC = matthews_corrcoef(val_labels, predictions)
        # #ROC = roc_auc_score(val_labels, predict_proba)
        # scores = [sn, sp, acc, MCC]
        # all_scores[i] = scores
    # print(all_scores)
    # with open('Rse18.txt', 'w') as file:
    #     # 使用字符串表示写入
    #     file.write(str(all_scores))
# ___________________________________训练、测试、独立________________________________________________________________
#
    #criterion = nn.CrossEntropyLoss()
    #criterion = nn.BCELoss()
   #  criterion = nn.BCEWithLogitsLoss()  # 二分类交叉熵损失
   #  optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
   #
   #  scheduler = ReduceLROnPlateau(optimizer, mode='min', factor=0.95, patience=2, verbose=True)
   #
   # # device = torch.device("cuda: 1" if torch.cuda.is_available() else "cpu")
   # # device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
   #  torch.cuda.set_device(1)
   #  device = torch.device("cuda")
   #  model.to(device)
   #  num_epochs = 40
   #
   #  for epoch in range(num_epochs):
   #      model.train()
   #      train_loss = 0
   #      correct = 0
   #      total = 0
   #
   #      for data, labels in train_loader:
   #          data, labels = data.to(device), labels.to(device)
   #          optimizer.zero_grad()
   #          outputs = model(data)
   #          loss = criterion(outputs, labels)
   #          loss.backward()
   #          optimizer.step()
   #
   #          train_loss += loss.item()
   #         # _, predicted = outputs.max(1)
   #          predicted = (outputs > 0.5).float()
   #          total += labels.size(0)
   #          correct += predicted.eq(labels).sum().item()
   #
   #      val_acc,val_loss,val_predict = evaluate(model, val_loader,device)
   #      test_acc,test_loss,test_predict = evaluate(model, test_loader,device)
   #      save_best_model(model, val_acc, test_acc)
   #      scheduler.step(val_loss)
   #      print(
   #          f"Epoch [{epoch + 1}/{num_epochs}], Loss: {train_loss / len(train_loader):.4f}, Accuracy: {100. * correct / total:.2f}%")
   #      print(f"验证集准确率: {val_acc:.2f}%, Loss:{val_loss:.4f}")
   #      print(f"测试集准确率: {test_acc:.2f}%, Loss:{test_loss:.4f}")

