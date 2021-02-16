---
title: "Installing Conda in Google Colab"
date: 2021-02-16
lastmod: 2021-02-16
draft: false
tags: ["Python", "Bioinformatics", "Colab", "Linux"]
categories: ["Python", "Bioinformatics", "Colab", "Linux"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

**Conda** is the environment and package management solution for a number of popular data science tools including **Pandas**, **Scikit-Learn**, **PyTorch**, **NVIDIA Rapids** and many others. 

**Conda** also dramatically simplifies the process of installing popular deep learning tools like **TensorFlow**.


Then, how to install Conda when using Google Colab?

<!--more-->

### Preliminaries

- To confirm which Python is being used by default in Google Colab
  
  ```python
    !which python
    #Output:
    #/usr/local/bin/python
  ```

- Check the version number of this default Python.

    ```python
    !python --version
    #Output:
    #Python 3.6.9
    ```

- Check to see if the PYTHONPATH variable has been set.

    ```python
    !echo $PYTHONPATH
    #Output:
    #/env/python
    ```

- Unset the PYTHONPATH variable

    ```python
    %env PYTHONPATH=
    #Output:
    #env: PYTHONPATH=
    ```

### Installing Miniconda

```python
%%bash
MINICONDA_INSTALLER_SCRIPT=Miniconda3-4.5.4-Linux-x86_64.sh
MINICONDA_PREFIX=/usr/local
wget https://repo.continuum.io/miniconda/$MINICONDA_INSTALLER_SCRIPT
chmod +x $MINICONDA_INSTALLER_SCRIPT
./$MINICONDA_INSTALLER_SCRIPT -b -f -p $MINICONDA_PREFIX
```

- Check Conda install and version

    ```python
    !which conda
    !conda --version
    #Output:
    #/usr/local/bin/conda
    #conda 4.5.4
    ```

- Miniconda has actually installed a slightly different version of Python.

    ```python
    !python --version
    #Output:
    #Python 3.6.5 :: Anaconda, Inc.
    ```

### Updating Conda

- Update Conda and all its dependencies to their most recent versions without updating Python to 3.7

    ```python
    %%bash
    conda install --channel defaults conda python=3.6 --yes
    conda update --channel defaults --all --yes
    ```

- New version of conda
  
    ```python
    !conda --version
    #Output:
    #conda 4.9.2
    ```

- Also, note that the Python version has changed yet again.

    ```python
    !python --version
    #Output:
    #Python 3.6.12 :: Anaconda, Inc.
    ```

### Appending to the sys.path

- Now that you have installed Miniconda you need to add the directory where Conda will install packages to the list of directories that Python will search when looking for modules to import. You can see the current list of directories that Python will search when looking for modules to import by inspecting the `sys.path`.

    ```python
    import sys
    sys.path
    ```

- The preinstalled packages included with Google Colab are installed into the `/usr/local/lib/python3.6/dist-packages` directory.

- To get an idea of what packages are available by simply listing the contents of this directory with Google Colab are installed

    ```python
    !ls /usr/local/lib/python3.6/dist-packages
    ```

- Any package that you install with Conda will be installed into the directory `/usr/local/lib/python3.6/site-packages`

- So you will need to add this directory to `sys.path` in order for these packages to be available for import.

    ```python
    import sys
    _ = (sys.path
            .append("/usr/local/lib/python3.6/site-packages"))

    ```

### Installing packages

- Include the `--yes` flag when installing your packages to avoid getting prompted to confirm the package plan.

    ```python
    !conda install --channel conda-forge featuretools --yes
    ```

### Summary

- In this article, I covered the process that I use to install and configure Miniconda when I need to use Conda to manage packages on Google Colab.
