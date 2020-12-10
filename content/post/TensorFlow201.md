---
title: "TensorFlow2 class1: Iris"
date: 2020-12-02
lastmod: 2020-12-10
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

TensorFlow2 notebook: Class1 神经网络计算.

<!--more-->

## In summary

### Prerequisites
- `TensorFlow2.3.0` + `Python3.6.9`

### 人工智能三学派
- 行为主义
- 符号主义
- 连接主义
  - 神经网络
  - 计算机仿真神经网络连接关系
    - 准备数据:'特征-标签'数据
    - 搭建网络:神经网络结构
    - 优化参数:训练网络获取最佳参数
    - 应用网络:输出分类，或输出预测结果

### 神经网络设计过程
- 分类目标:
  - 0狗尾草Iris
  - 1杂色Iris: 花萼长>花萼宽 and 花瓣长/花瓣宽>2
  - 2弗吉尼亚Iris
- 采集大量数据对作为数据集: 输入特征(花萼长，花萼宽，花瓣长，花瓣宽), 人工标定的标签(对应的类别)
- 将已有数据集喂入搭建好的神经网络结构
- 随机初始化所有参数, 然后反向传播进行网络优化参数,得到模型
  - 损失函数输出最小(预测值与标准答案之间的差距)，得到所有的最优参数
    - 梯度下降法
      - 设置合适的学习率(超参数): 不可以过大或者过小
- 读入新输入特征(待测)
- 输出识别结果: 所属分类

### 张量生成
  
### 常用函数

- reduce_mean,reduce_sum
  - axis=0:纵向，经度方向
  - axis=1:横向，维度方向
- 维度相同的张量可以做四则运算
  - add, subtract, multiply, divide.
- 平方, 次方, 开方
  - square, pow, sqrt
- 矩阵乘
  - tf.matmul
- 将(特征, 标签)进行配对
  
  ```python
  tf.data.Dataset.from_tensor_slices((features, labels))
  ```

- 梯度
  
  ```python
  tf.GradientTape
  ```

- 枚举, python内置函数
  
  ```python
  enumerate
  ```

- 独热编码,作为标签
  - tf.one_hot
  - 1: 表示是
  - 0: 表示非
- 使输出符合概率分布 -tf.nn.softmax -->n分类n输出调用softmax()函数
- 自减
  
  ```python
  assign_sub
  ```

- 张量方向最大值的索引
  
  ```python
  tf.argmax()
  ```

### Iris 数据集读入

### 神经网络实现iris分类

Attach is the file of [TensorFlow2class1.ipynb](TensorFlow2class1.ipynb), or view it via the [link](https://colab.research.google.com/drive/1IslgEamApGpdgwBK2YW0zgGBJJ9L7VPb?usp=sharing).

Show me the code <i class="far fa-hand-point-down"></i>

  [![fig1](fig1.png)](https://gist.github.com/xiaonilee/92768a4f7c91fc7553e58f1b2a7e02d5)
  [![fig2](fig2.png)](https://gist.github.com/xiaonilee/92768a4f7c91fc7553e58f1b2a7e02d5)

{{< gist 92768a4f7c91fc7553e58f1b2a7e02d5 >}}


  