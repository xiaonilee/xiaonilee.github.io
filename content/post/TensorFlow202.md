---
title: "TensorFlow2 class2: 数据网络优化"
date: 2020-12-03
lastmod: 2021-02-19
draft: false
tags: ["TensorFlow2", "Machine Learning"]
categories: ["TensorFlow2", "Machine Learning", "big data"]
author: "Xiaoni"

weight: 1

mathjax: true

# menu:
#   main:
#     parent: "docs"
#     weight: 1
---

TensorFlow2 notebook: Class2 神经网络的优化方法, 学习率,激活函数,损失函数和正则化的使用, 用Python语言写出5种反向传播优化器(SGD, Momentum, Adagrad, RMSProp, Adam).

<!--more-->

## In summary

### 2.1 预备知识
- `tf.where(条件语句, 真返回A, 假返回B)`
- `np.random.RandomState.rand()`, 返回一个[0,1)之间的随机数
- `np.vstack()`, 将两个数组按照垂直方向叠加
- `np.mgrid[起始值:结束值:步长, 起始值:结束值:步长]; [起始值:结束值)`
- `.ravel()`, 将x变为一维数组，把.前变量拉直
- `np.c_[]`, 使得返回的间隔数值点配对  

### 2.2 神经网络复杂度  
- 神经网络(NN)复杂度:多用NN层数和NN参数的个数表示
  - 空间复杂度, 层数 = 隐藏层的层数 + 1个输出层
  - 时间复杂度, 乘加运算次数
- 指数衰减学习率
  - 初始学习率 * 学习率衰减率 ^ (当前轮数 / 多少轮衰减一次)

### 2.3 激活函数
- 非线性函数
- 大大提升列模型的表达力
- 优秀的激活函数特征
  - 非线性, 多层神经网络可逼近所有函数
  - 可微性, 优化器大多用梯度下降更新参数
  - 单调性, 保证单层网络的损失函数是凸函数
  - 近似恒等性: f(x)≈x, 神经网络更稳定
- Sigmoid函数

  ```python
  tf.nn. sigmoid(x)
  Tanh函数
  tf.math. tanh(x)
  Relu函数,对于初学者首选该函数
  tf.nn.relu(x)
  Leaky Relu函数
  tf.nn.leaky_relu(x)
  ```


### 2.4 损失函数
- 预测值(y)与已知答案(y_)的差距
- NN 优化的终极, 重点目标, 使得loss最小
  - MSE

  ```python
  loss_mse = tf.reduce_mean()
  ```

  - 自定义
    
    ```python
    tf.reduce_sum(tf.where(tf.greater())
    ```

  - CE
    
    ```python
    tf.losses.categorical_crossentropy()
    ```

- 案例
  - 预测酸奶日销量y，x1、x2是影响日销量的因素
  - 没有真实数据集，现通过随机生成来拟造 
  
### 2.5 缓解**过拟合**
- 欠拟合, 学习不够彻底
  - 增加输入特征项
  - 增加网络参数
  - 减少正则化参数
- 过拟合, 缺乏泛化力
  - 数据清洗，减少集中噪声
  - 增大训练集
  - 采用正则化
  - 增大正则化参数

### 2.6 优化器更新网络参数 
- 优化器
  - ***引导神经网络更新参数的工具***
- 待优化参数𝒘, 损失函数loss, 学习率𝒍r, 每次迭代一个batch(2^n), t表示当前batch迭代的总次数.
  - 一阶动量:与梯度相关的函数
  - 二阶动量:与梯度平方相关的函数
- 优化器
  - SGD(无momentum), 常用的梯度下降法
  - SGDM(含momentum的SGD), 在SGD基础上增加一阶动量
  - Adagrad，在SGD基础上增加二阶动量
  - RMSProp，SGD基础上增加二阶动量
  - Adam, 同时结合SGDM一阶动量和RMSProp二阶动量

Attach is the file of [TensorFlow2class2.ipynb](TensorFlow2class2.ipynb), or view it via the [link](https://colab.research.google.com/drive/19S0UpKcWc_l6Xi-aNcu-Tr3OaNi7SyTX?usp=sharing).

Show me the code <i class="far fa-hand-point-down"></i>

[![fig1](fig1.png)](https://gist.github.com/xiaonilee/d2854ec53cfb4c8797006393af95c1f8)

{{< gist d2854ec53cfb4c8797006393af95c1f8 >}}
