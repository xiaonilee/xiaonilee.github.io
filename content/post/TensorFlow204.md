---
title: "TensorFlow2 class4: 网络八股扩展"
date: 2020-12-08
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

TensorFlow2 notebook: Class4: 网络八股进行扩展, 增加自制数据集, 数据增强, 断点续训, 参数提取和acc/loss可视化, 最后实现给图识物的应用程序.

<!--more-->

## In summary

### 搭建网络八股总览
- 自制数据集，解决本领域应用
  - 将自己本领域的数据和标签赋值给`x_train`, `y_train`, `x_test`, `y_test`
- 数据增强，扩充数据集
  - 数据增强的代码: 数据量过少
- 断点续训，存取模型
  - 断点续训的代码, 实时保存最优模型, 从而无需总是从零开始训练模型
- 参数提取的代码, 将参数存入文本
- acc/loss曲线绘制的代码, 可视化训练效果
- 应用程序, 给图识物

### 自制数据集
- 替代**TensorFlow**中自带函数`.load_data()`

### 数据增强
- TensorFlow中的数据增强函数

```python
tf.keras.preprocessing.image.ImageDataGenerator()
```

- 调整函数中的参数
  - `rescale=` 所有数据将乘以该数值
  - `rotation_range=`  随机旋转角度数范围
  - `width_shift_range=` 随机宽度偏移量
  - `height_shift_range=` 随机高度偏移量
  - `horizontal_flip=` 是否随机水平翻转
  - `zoom_range=` 随机缩放的范围[1-n，1+n]
  
### 断点续训
- 读取模型
  - TensorFlow中函数`model.load_weights(路径文件名)`,读取已有模型参数
- 保存模型

    ```python
    cp_callback=tf.keras.callbacks.ModelCheckpoint()
    ```

  - 结果保存在checkpoint文件夹中
    - [checkpoint](checkpoint)
    - [mnist.ckpt.data-00000-of-00001](mnist.ckpt.data-00000-of-00001)
    - [mnist.ckpt.index](mnist.ckpt.index)
- 继续从上次模型训练结果进行
- 生成的`checkpoint`文件夹中存放的就是模型参数

### 参数提取
- 提取可训练参数
  - model.trainable_variables 返回模型中可训练的参数
- 设置print输出格式

    ```python
    np.set_printoptions(threshold=超过多少省略显示)
    threshold=np.inf表示无限大
    ```

- for loop将所有可训练的参数存入到[weights.txt](weights.txt)文件中

### acc&loss可视化
- history:
  - 训练集loss: loss
  - 测试集loss: val_loss
  - 训练集准确率: sparse_categorical_accuracy
  - 测试集准确率: val_sparse_categorical_accuracy
- history.history从model.fit中提取

    ```python
    acc = history.history['sparse_categorical_accuracy']
    val_acc = history.history['val_sparse_categorical_accuracy']
    loss = history.history['loss']
    val_loss = history.history['val_loss']
    ```

### 给图识物
- 复现模型

    ```python
    model = tf.keras.models.Sequential([
    tf.keras.layers.Flatten(),
    tf.keras.layers.Dense(128, activation='relu'), 
    tf.keras.layers.Dense(10, activation='softmax’)])
    ```

- 加载保存在checkpoint文件中的参数

    ```python
    model.load_weights(model_save_path)`
    ```

- 预测结果
  
  ```python
  result = model.predict(x_predict)
  ```


Attach is the file of [TensorFlow2class4.ipynb](TensorFlow2class4.ipynb), or view it via the [link](https://colab.research.google.com/drive/1lbZj25hoCJzQbwSJRts0X1gwdLW-Diw7?usp=sharing).

Show me the code <i class="far fa-hand-point-down"></i>

[![fig](fig1.png)](https://gist.github.com/xiaonilee/65f56d831a260679125da6e54872af00)

{{< gist 65f56d831a260679125da6e54872af00 >}}
