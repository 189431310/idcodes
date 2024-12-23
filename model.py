import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import random


class BasicBlock(nn.Module):
    def __init__(self, inchannel, outchannel, stride=1):
        super(BasicBlock, self).__init__()
        self.left = nn.Sequential(
            nn.Conv2d(inchannel, outchannel, kernel_size=3, stride=stride, padding=1, bias=False),
            nn.BatchNorm2d(outchannel),
            nn.ReLU(inplace=True),
            nn.Conv2d(outchannel, outchannel, kernel_size=3, stride=1, padding=1, bias=False),
            nn.BatchNorm2d(outchannel)
        )
        self.shortcut = nn.Sequential()
        if stride != 1 or inchannel != outchannel:
            self.shortcut = nn.Sequential(
                nn.Conv2d(inchannel, outchannel, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(outchannel)
            )

    def forward(self, x):
        out = self.left(x)
        out += self.shortcut(x)
        out = F.relu(out)
        return out

class ResNet_18(nn.Module):
    def __init__(self, ResidualBlock, num_classes=2):  # 修改为二分类 num_classes=2
        super(ResNet_18, self).__init__()
        self.inchannel = 64
        self.conv1 = nn.Sequential(
            nn.Conv2d(1, 64, kernel_size=3, stride=1, padding=1, bias=False),  # 修改为单通道输入 (1代表输入通道数)
            nn.BatchNorm2d(64),
            nn.ReLU(),
        )
        self.layer1 = self.make_layer(ResidualBlock, 64, 2, stride=1)
        self.layer2 = self.make_layer(ResidualBlock, 128, 2, stride=2)
        self.layer3 = self.make_layer(ResidualBlock, 256, 2, stride=2)
        self.layer4 = self.make_layer(ResidualBlock, 512, 2, stride=2)
        self.fc = nn.Linear(512, num_classes)

    def make_layer(self, block, channels, num_blocks, stride):
        strides = [stride] + [1] * (num_blocks - 1)
        layers = []
        for stride in strides:
            layers.append(block(self.inchannel, channels, stride))
            self.inchannel = channels
        return nn.Sequential(*layers)

    def forward(self, x):
        out = self.conv1(x)
        out = self.layer1(out) #(64,64,101,72)
        out = self.layer2(out) #(64,128,51,36)
        out = self.layer3(out) #(64,128,26,18)
        out = self.layer4(out) #(64,512,13,9)
        out = F.avg_pool2d(out, (13,9))  # 萃取全局特征 (64,512,3,2)
        out = out.view(out.size(0), -1)
        out = self.fc(out)
        return out



def set_seed(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def evaluate(model, loader,device):
    model.eval()
    correct = 0
    total = 0
    loss = 0  # Accumulate loss
    #criterion =  nn.BCEWithLogitsLoss()   # Loss function
    criterion = nn.CrossEntropyLoss()
    #criterion = nn.BCELoss()

    result ={"predicted_label":[], "predicted_probability":[]}

    with torch.no_grad():
        for data, labels in loader:
            data, labels = data.to(device), labels.to(device)
            outputs = model(data)
            loss += criterion(outputs, labels).item()
            _,predicted = outputs.max(1)
            #predicted =(outputs>0.5).float()
            result['predicted_label'].extend(predicted.cpu().numpy().tolist())
            #result['predicted_probability'].extend(outputs.cpu().numpy().tolist())
            total += labels.size(0)
            correct += predicted.eq(labels).sum().item()
    accuracy = 100. * correct / total
    avg_loss = loss / len(loader)

    return accuracy, avg_loss, result


def set_seed(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


best_model_info = {'val_acc': 0, 'test_acc': 0, 'model_state': None}

# 保存模型函数
def save_best_model(model, val_acc, test_acc, save_path='CNN_best_model.pth'):
    global best_model_info
    if val_acc > 80 and test_acc > 80 and (val_acc > best_model_info['val_acc'] or test_acc > best_model_info['test_acc']):
        print(f"保存新模型: val_acc={val_acc:.2f}, test_acc={test_acc:.2f}")
        best_model_info['val_acc'] = val_acc
        best_model_info['test_acc'] = test_acc
        best_model_info['model_state'] = model.state_dict()  # 保存模型权重
        torch.save(model.state_dict(), save_path)