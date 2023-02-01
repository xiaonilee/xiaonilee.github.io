---
title: "win10系统下Rstudio中的乱码问题"
date: 2023-01-25
lastmod: 2023-01-25
draft: false
tags: ["program", "R", "win10"]
categories: ["program", "win10"]
author: "Xiaoni"

weight: 1

mathjax: true

# menu:
#   main:
#     parent: "docs"
#     weight: 1
---

既往使用mac系统，未曾见过乱码现象。然而，目前工作（同事）主要使用win10系统。

对于来自不同同事的R脚本，会出现各种的乱码。期间尝试过各种解决办法，但均没有从根本上解决。

在此记录这个好的解决办法，再也不见同类问题。


<!--more-->

### 问题
- 如图1:
  - ![code1](code1.png)

### 解决办法
- 在开头运行Sys.setlocale('LC_CTYPE', 'Chinese' )。
- 重启 R。
- 加入代码：encoding = 'UTF-8'。
- 如图2：
  - ![code2](code2.png)
- Done.