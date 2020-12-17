---
title: "TensorFlow2 class6: 循环神经网络"
date: 2020-12-13
lastmod: 2020-12-17
draft: false
tags: ["TensorFlow2", "ML"]
categories: ["TensorFlow2", "ML", "big data"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

TensorFlow2 notebook: Class6: 循环神经网络.

<!--more-->

## In summary

### Prerequisites
- Setup virtul environment
- `python3.6.9` + `TensorFlow2.3.0` + `sklearn` + `pandas` + `matplotlib`

### 6.0 回顾
- **卷积神经网络:**
  - 卷积
    - 特征提取器, CBAPD
  - 卷积神经网络
    - 借助**卷积核**提取**空间特征**后，送入全连接网络, 实现离散数据的分类

### 6.1 循环核
- 具有记忆力
- **不同时刻的参数**共享, 实现了对时间序列的信息提取
- $$输入指定维度的x_t\ 输出指定维度的y_t\ $$
- $$n个记忆体内存储的状态信息h_t\  + 下面待训练参数矩阵 w_{xh}\  + 周围侧面待训练参数矩阵 w_{hh}\ + 上面待训练参数矩阵w_{hy}\ $$
  - 前向传播时: $$h_t\ 在每个时刻都被刷新，w_{xh}\ w_{hh}\ w_{hy}\ 自始至终都固定不变$$
  - 反向传播时: $$x_h\ w_{hh}\ w_{hy}\ 被梯度下降法更新$$
  - $$h_t = 激活函数tanh(x_tw_{xh} + h_{t-1}w_{hh} + bh)\ $$
  - $$y_t = 激活函数softmax(h_tw_{hy}+by)\ 是循环网络末层, 是一层全连接$$

### 6.2 循环核时间步展开
- 将循环核按时间步展开, 按照时间轴方向展开t步: $$x_0\ \ x_1\ \ x_2\ \ ... x_{t-1} \ \ x_t\ $$
- $$每个时刻h_t 被刷新,记忆体周围的参数矩阵w_{xh}  + w_{hh} + w_{hy} 固定不变(类似于人脑中的记忆体, 每个时刻会根据当前的输入而更新, 当前时刻的推理是根据以往的知识积累, 用固化下来的参数矩阵进行推理判断), 进行训练优化, 得到最好的参数矩阵, 进行前向传播, 输出连续数据的预测结果\$$

### 6.3 循环计算层
- **每个**循环核构成一**层**循环计算层
  - 每个循环核中**记忆体的个数**, 可以根据需求**任意指定**
- 层数是向输出方向增长的

### 6.4 TF描述循环计算层

```python
tf.keras.layers.SimpleRNN(记忆体个数，activation=‘激活函数’, 
return_sequences=是否每个时刻输出ht到下一层)
activation=‘激活函数’ (不写，默认使用tanh) 
return_sequences=True 循环核各时刻会把ht推送到到下一层
return_sequences=False 循环核仅在最后一个时刻把ht推送到到下一层(默认)

例:SimpleRNN(3, return_sequences=True)
```

- 送入RNN时， 要求x_train维度是3维:
  - [送入样本个数, 循环核时间展开步的步数, 每个时间步输入特征个数]

### 6.5 循环计算过程I
- 字母预测
  - 目标: 输入a预测出b，输入b预测出c， 输入c预测出d，输入d预测出e，输入e预测出a
  - 将a b c d e用数字表示出来: 独热码进行编码
    - a: 10000
    - b: 01000
    - c: 00100
    - d: 00010
    - e: 00001
  - $$随机生成w_{xh}  \ \ w_{hh} \ \ w_{hy} 三参数的矩阵值$$
  - 记忆体个数设为: $$3$$
  - 上一时刻, 起始时刻的记忆体信息: $$ h_{t-1}=0$$
  - $$代入随机生成的三参数数值和初始h_{t-1}值,计算当前时刻的h_t = tanh(x_tw_{xh} + h_{t-1}w_{hh} + bh), 将记忆体存储的状态信息进行更新$$
  - $$将得到的h_t矩阵值代入公式softmax(h_tw_{hy}+by) \\ 预测出y_t结果, 概率最大的即为结果\$$

### 6.6 字母预测onehot_1pre1
- RNN实现连续输入1个字母预测下一个字母
- 字母进行独热码编码one hot
- [rnn_onehot_1pre1.ckpt.index](606/rnn_onehot_1pre1.ckpt.index)
- [rnn_onehot_1pre1.ckpt.data-00000-of-00001](606/rnn_onehot_1pre1.ckpt.data-00000-of-00001)
- [checkpoint](606/checkpoint)
- [weights](606/weights.txt)

### 6.7 循环计算过程II
- 连续输入4个(多个)字母, 预测下一个字母
- 使用3个记忆体, 初始时刻为[0 0 0], 每个时间步之后被计算依次更新 
- 每个时间步具有相同的$$w_{xh}  \ \ w_{hh} \ \ bh $$
  
```r
# 矩阵计算过程
xtb <- matrix(c(0, 1, 0, 0, 0), nrow=1, ncol=5, byrow=TRUE)
xtc <- matrix(c(0, 0, 1, 0, 0), nrow=1, ncol=5, byrow=TRUE)
xtd <- matrix(c(0, 0, 0, 1, 0), nrow=1, ncol=5, byrow=TRUE)
xte <- matrix(c(0, 0, 0, 0, 1), nrow=1, ncol=5, byrow=TRUE)

ht1 <- matrix(c(-0.9, 0.2, 0.2), nrow=1, ncol = 3, byrow = TRUE)
ht2 <- matrix(c(0.8, 1.0, 0.8), nrow=1, ncol = 3, byrow = TRUE)
ht3 <- matrix(c(0.6, 0.5, -1.0), nrow=1, ncol = 3, byrow = TRUE)
ht4 <- matrix(c(-1.0, -1.0, 0.8), nrow=1, ncol = 3, byrow = TRUE)

whh <- matrix(c(-0.9,-0.9,-0.9,0.5,0.9,-0.3,1.0,0.3,-1.5), 
              nrow = 3, ncol = 3, byrow = TRUE)

wxh <- matrix(c(1.2,-1.3,1.1,-1.5,0.2,0.3,-0.3,1.7,0.7,-0.1,0.1,-0.1,-1.2,-1.5,0.3), 
              nrow = 5, ncol = 3, byrow = TRUE)

wyt <- matrix(c(-1.3,0.5,-0.7,-0.2,0.8, -1.4,-0.8,-1.2,0.9,1.4,0.7,1.1,-1.2,1.3,-1.1), 
              nrow = 3, ncol = 5, byrow = TRUE)

xtb %*% wxh
xtc %*% wxh
xtd %*% wxh
xte %*% wxh

round(ht1 %*% whh, 1)
round(ht2 %*% whh, 1)
round(ht3 %*% whh, 1)
round(ht4 %*% wyt, 1)
```

### 6.8 字母预测onehot_4pre1
- RNN实现连续输入4个字母预测下一个字母
  - 输入abcd输出e 
  - 输入bcde输出a 
  - 输入cdea输出b 
  - 输入deab输出c 
  - 输入eabc输出d
- 字母进行独热码编码one hot

#### 6.8.1 交互式输出字母预测结果
- [rnn_onehot_4pre1.ckpt.index](608_1/rnn_onehot_4pre1.ckpt.index)
- [rnn_onehot_4pre1.ckpt.data-00000-of-00001](608_1/rnn_onehot_4pre1.ckpt.data-00000-of-00001)
- [checkpoint](608_1/checkpoint)
- [weights](608_1/weights.txt)

#### 6.8.2 一次性输出字母预测结果
- [rnn_onehot_4pre1.ckpt.index](608_2/rnn_onehot_4pre1.ckpt.index)
- [rnn_onehot_4pre1.ckpt.data-00000-of-00001](608_2/rnn_onehot_4pre1.ckpt.data-00000-of-00001)
- [checkpoint](608_2/checkpoint)
- [weights](608_2/weights.txt)

### 6.9 Embedding编码
- 前面用到的独热码位宽要与词汇量一致
- Embeddig是一种单词编码方法
- 用低维向量实现了编码
- 编码通过神经网络训练优化，能表达出单词间的相关性
- tf函数表示

  ```python
    tf.keras.layers.Embedding(词汇表大小，编码维度)
  ```

- x_train 要求维度格式

  ```python
  [送入样本数， 循环核时间展开步数]
  ```

### 6.10 字母预测Embedding_1pre1
- RNN实现连续输入1个字母预测下一个字母
- Embedding编码

#### 6.10.1 交互式输出字母预测结果
- [run_embedding_1pre1.ckpt.index](610_1/run_embedding_1pre1.ckpt.index)
- [run_embedding_1pre1.ckpt.data-00000-of-00001](610_1/run_embedding_1pre1.ckpt.data-00000-of-00001)
- [checkpoint](610_1/checkpoint)
- [weights](610_1/weights.txt)

#### 6.10.2 一次性输出字母预测结果
- [run_embedding_1pre1.ckpt.index](610_2/run_embedding_1pre1.ckpt.index)
- [run_embedding_1pre1.ckpt.data-00000-of-00001](610_2/run_embedding_1pre1.ckpt.data-00000-of-00001)
- [checkpoint](610_2/checkpoint)
- [weights](610_2/weights.txt)

### 6.11 字母预测Embedding_4pre1
- RNN实现连续输入4个字母预测下一个字母
- Embedding编码
- 词汇量26个
- 字母a->z用数字0-25依次一一进行表示
- for loop 构建`x_train` `y_train`
  - 从数字列表中把连续4个数作为输入特征x_train
  - 第5个数作为标签y_train

#### 6.11.1 交互式输出字母预测结果
- [rnn_embedding_4pre1.ckpt.index](611_1/rnn_embedding_4pre1.ckpt.index)
- [rnn_embedding_4pre1.ckpt.data-00000-of-00001](611_1/run_embedding_4pre1.ckpt.data-00000-of-00001)
- [checkpoint](611_1/checkpoint)
- [weights](611_1/weights.txt)

#### 6.11.2 一次性输出字母预测结果
- [rnn_embedding_4pre1.ckpt.index](611_2/rnn_embedding_4pre1.ckpt.index)
- [rnn_embedding_4pre1.ckpt.data-00000-of-00001](611_2/run_embedding_4pre1.ckpt.data-00000-of-00001)
- [checkpoint](611_2/checkpoint)
- [weights](611_2/weights.txt)

### 6.12 RNN实现股票预测
- 数据源[SH600519.csv](SH600519.csv)
- 数据下载源码p37_tushare.py
  
  ```python
  import tushare as ts
  import matplotlib.pyplot as plt

  df1 = ts.get_k_data('600519', ktype='D', start='2010-04-26', end='2020-04-26')

  datapath1 = "./SH600519.csv"
  df1.to_csv(datapath1)
  ```

- 结果文件
  - [rnn_stock.ckpt.index](612/rnn_stock.ckpt.index)
  - [rnn_stock.ckpt.data-00000-of-00001](612/rnn_stock.ckpt.data-00000-of-00001)
  - [checkpoint](612/checkpoint)
  - [weights](612/weights.txt)
### 6.00 半讲小结
- 上述传统循环网络RNN
  - 通过记忆体实现短期记忆, 来进行连续数据的预测 
  - 若连续数据序列变长时, 会使**展开时间步**过长, 在反向传播更新参数时,梯度要按照**时间步连续相乘**而导致梯度消失

### 6.13 LSTM实现股票预测(LSTM计算过程_TF描述LSTM层)
- 于1997年提出，通过门控单元改善了RNN长期依赖问题
- 长短记忆网络LONG SHORT-TERM MEMORY
- **3门限**
  - 输入门(门限)
  - 遗忘门(门限)
  - 输出门(门限)
  - $$当前时刻输入特征x_t, 上一时刻的短期记忆函数h_{t-1},输入门 i_t, 遗忘门f_t, 输出门o_t三门中各自待训练参数矩阵W_i, W_f, w_o, 与相应的待训练偏置项b_i,b_f,b_o,\\ 各自经过sigmoid函数$$
  - 门限范围0-1
- **细胞态**,由过去的长期记忆+现在的候选态构成
  - $$表示长期记忆c_t$$
  - 上一时刻的长期记忆(过去的长期记忆) x 遗忘门 + 当前时刻归纳出的新知识(候选态) x 输入门
- **记忆体**的输出, 是长期记忆的一部分
  - $$表示短期记忆h_t$$
  - 细胞态经过tanh激活函数后(留存在脑中的长期记忆) x 输出门
- **候选态**
  - 等待存入长期记忆, 归纳出的新知识
- TF描述LSTM层

  ```python
  tf.keras.layers.LSTM(记忆体个数，return_sequences=是否返回输出)
  return_sequences=True 各时间步输出ht,多用于中间层网络
  return_sequences=False 仅最后时间步输出ht(默认), 多用于最后一层网络
  ```

- 结果文件
  - [LSTM_stock.ckpt.index](613/LSTM_stock.ckpt.index)
  - [LSTM_stock.ckpt.data-00000-of-00001](613/LSTM_stock.ckpt.data-00000-of-00001)
  - [checkpoint](613/checkpoint)
  - [weights](613/weights.txt)
### 6.14 GRU实现股票预测(GRU计算过程_TF描述GRU层)
- GRU网络于2014年提出，优化LSTM结构
- 两门
  - $$更新门z_t$$
  - $$重置门r_t$$
  - 门限范围0-1
- $$ 记忆体h_t融合了长期记忆和短期记忆$$
- $$ 上一时刻记忆体 x (1-更新门) + 候选隐藏层 x 更新门$$
- $$候选隐藏层取决于上一时刻记忆体 * 重置门以及当前输入特征x_t$$

- TF描述GRU层

  ```python
  tf.keras.layers.GRU(记忆体个数，return_sequences=是否返回输出) return_sequences=True 各时间步输出ht
  return_sequences=False 仅最后时间步输出ht(默认)
  ```

- 结果文件
  - [stock.ckpt.index](614/stock.ckpt.index)
  - [stock.ckpt.data-00000-of-00001](614/stock.ckpt.data-00000-of-00001)
  - [checkpoint](614/checkpoint)
  - [weights](614/weights.txt)

Attach is the file of [TensorFlow2class6.ipynb](TensorFlow2class6.ipynb), or view it via the [link](https://colab.research.google.com/drive/1NJgM3NG0igL-zTx3k-EkNk0-7rR5ES9K?usp=sharing).

Show me the code <i class="far fa-hand-point-down"></i>

[![fig2](fig2.png)](https://gist.github.com/xiaonilee/0f5034da455136055559511a9bfe23e9)
[![fig3](fig3.png)](https://gist.github.com/xiaonilee/0f5034da455136055559511a9bfe23e9)

{{< gist 0f5034da455136055559511a9bfe23e9 >}}
