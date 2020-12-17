---
title: "TensorFlow2 class3: 神经网络八股"
date: 2020-12-07
lastmod: 2020-12-10
draft: false
tags: ["TensorFlow2", "ML"]
categories: ["TensorFlow2", "ML", "big data"]
author: "Xiaoni"

weight: 1

mathjax: true

# menu:
#   main:
#     parent: "docs"
#     weight: 1
---

TensorFlow2 notebook: Class3 神经网络八股.

<!--more-->

## In summary

### 搭建神经网络八股sequential
- **6步法**搭建神经网络(也是提纲):
  - 3.1.1 用TensorFlow API: tf.keras搭建网络八股
    - import
      - 所需模块
  - 3.1.2 指定train和test的特征以及标签(也可以将test数据集在第5步fit中指定)
  - 3.1.3 model = tf.keras.models.Sequential
    - Sequential([**一个封装从输入层到输出层网络结构的容器**])中执行
    - 搭建网络结构, 逐层网络进行描述，equal前向传播
  - 3.1.4 model.compile，配置神经网络的训练方法
    - compile()中执行
    - 所选optimizer优化器
    - 所选loss函数
    - 所选metrics评测指标, acc
  - 3.1.5 model.fit()
    - fit()中执行
    - train, test的输入特征和标签
    - 选定每次喂入神经网络的batch_size数目
    - 选定epochs数目
    - validation_split
    - validation_freq: 多少次epoch测试集验证一次准确率
  - 3.1.6 model.summary
    - summary()中执行
    - 打印网络结构和参数统计

### 搭建神经网络八股class
- 带有跳连的非顺序网络结构
  - `class MyModel(Model) model=MyModel` 进行封装
  - `class MyModel(Model)`
    - 括号中的Model继承了TensorFlow中的Model类
    - `def init(self)`: super(MyModel, self).init(),定义网络结构块
    - `def call(self,x)`: 调用网络结构块，实现前向传播 return y
  - `model = MyModel()`

### MNIST数据集
- 像素点灰度值数据
- 由7万张(28x28) 像素点的0~9手写数字图片和标签构成, 其中
  - 6万张作为训练集
  - 1万张作为测试集
- 导入数据集

```python
mnist = tf.keras.datasets.mnist
(x_train, y_train) , (x_test, y_test) = mnist.load_data()
```

- 将数据拉伸为一维数组(共计**28x28=784**个数值), 作为输入特征，输入神经网络

### FASHION数据集
- 像素点灰度值数据
- 由7万张(28x28) 像素点的衣裤图片和标签构成, 其中
  - 6万张作为训练集
  - 1万张作为测试集
- 标签Label, 有10类:
    0 T-shirt/top
    1 Trouser
    2 Pullover
    3 Dress
    4 Coat
    5 Sandal
    6 Shirt
    7 Sneaker
    8 Bag
    9 Ankle boot
- 导入数据集

```python
fashion = tf.keras.datasets.fashion_mnist
(x_train, y_train) , (x_test, y_test) = fashion.load_data()
```

- 将数据拉伸为一维数组(共计28x28=784个数值), 作为输入特征，输入神经网络

Attach is the file of [TensorFlow2class3.ipynb](TensorFlow2class3.ipynb), or view it via the [link](https://colab.research.google.com/drive/12MY9fU_I2cXfFsq9mJ0Td0RdnooKLDg2?usp=sharing).

Show me the code <i class="far fa-hand-point-down"></i>

[![fig4](fig4.png)](https://gist.github.com/xiaonilee/2c2cc1f3da7560dd605184f24ea5c3bb)

[![fig5](fig5.png)](https://gist.github.com/xiaonilee/2c2cc1f3da7560dd605184f24ea5c3bb)

{{< gist 2c2cc1f3da7560dd605184f24ea5c3bb >}}
