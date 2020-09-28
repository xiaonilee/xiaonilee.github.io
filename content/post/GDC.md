---
title: "GDC Notebook"
date: 2020-09-25
lastmod: 2020-09-27
draft: false
tags: ["TCGA", "Bioinformatics", "GDC"]
categories: ["Cancer Research", "Database"]
author: "Xiaoni"

weight: 1

# You can also close(false) or open(true) something for this content.
# P.S. comment can only be closed
# comment: false
# toc: false

# You can also define another contentCopyright. e.g. contentCopyright: "This is another copyright."
contentCopyright: ''
# reward: false
mathjax: true
menu:
  main:
    parent: "docs"
    weight: 1
---

The Genomic Data Commons ([GDC](https://portal.gdc.cancer.gov/)), is a research program of the National Cancer Institute (NCI).

The **mission** of the GDC is to provide the cancer research community with a **unified** data repository that enables data sharing across cancer genomic studies in support of **precision medicine**.

<!--more-->

## In Brief

- Download Data via GDC
  
- Data Wrangling
  
## Instruction of Download

- Click on `Repository` and choose `Cases` to setup interested data: TCGA-LUNG.
  
  ![caseset](caseset.png)

- Choose `Files` and select interested demand, and then click on `Manifest` to download data1.
  
  ![filesset1](filesset1.png)

- Continue, select `Isoform Expression Quantification` to change **Data Type** and click on `Manifest` to download data2.
  
  ![filesset2](filesset2.png)

- Keep `Cases` selection, remove all **previously** selection on the page of `Files`.

- Select clinical from `Data Category` and **bcr xml** from `Data Format`, then Download data3.
  
  ![filesset3](filesset3.png)

- Count lines in downloaded data.

```markdown

$ wc -l gdc_manifest.2020-09-27-*

Output
     568 gdc_manifest.2020-09-27-LUNG-miRNA-isoform.txt
     568 gdc_manifest.2020-09-27-LUNG-miRNA-seq.txt
     523 gdc_manifest.2020-09-27-clinical.txt
    1659 total
```

- Download the [GDC Data Transfer Tool Client](https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.0_OSX_x64_1.zip)

```markdown

wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.0_OSX_x64_1.zip

unzip gdc-client_v1.6.0_OSX_x64_1.zip
```

- Common Usage of gdc-client

```markdown

$ ./gdc-client --help

output

commands:
  {download,upload,settings}
                        for more information, specify -h after a command
    download            download data from the GDC
    upload              upload data to the GDC
```

- Usage of gdc-client download.

```markdown
$ ./gdc-client download --help

output
usage: gdc-client download [-h] [--debug] [--log-file LOG_FILE] [--color_off]
                           [-t TOKEN_FILE] [-d DIR] [-s server]
                           [--no-segment-md5sums] [--no-file-md5sum]
                           [-n N_PROCESSES]
                           [--http-chunk-size HTTP_CHUNK_SIZE]
                           [--save-interval SAVE_INTERVAL] [-k]
                           [--no-related-files] [--no-annotations]
                           [--no-auto-retry] [--retry-amount RETRY_AMOUNT]
                           [--wait-time WAIT_TIME] [--latest] [--config FILE]
                           [-m MANIFEST]
                           [file_id [file_id ...]]

positional arguments:
  file_id               The GDC UUID of the file(s) to download

optional arguments:
  -h, --help            show this help message and exit
  --debug               Enable debug logging. If a failure occurs, the program
                        will stop.
  --log-file LOG_FILE   Save logs to file. Amount logged affected by --debug
  --color_off           Disable colored output
  -t TOKEN_FILE, --token-file TOKEN_FILE
                        GDC API auth token file
  -d DIR, --dir DIR     Directory to download files to. Defaults to current
                        directory
  -s server, --server server
                        The TCP server address server[:port]
  --no-segment-md5sums  Do not calculate inbound segment md5sums and/or do not
                        verify md5sums on restart
  --no-file-md5sum      Do not verify file md5sum after download
  -n N_PROCESSES, --n-processes N_PROCESSES
                        Number of client connections.
  --http-chunk-size HTTP_CHUNK_SIZE, -c HTTP_CHUNK_SIZE
                        Size in bytes of standard HTTP block size.
  --save-interval SAVE_INTERVAL
                        The number of chunks after which to flush state file.
                        A lower save interval will result in more frequent
                        printout but lower performance.
  -k, --no-verify       Perform insecure SSL connection and transfer
  --no-related-files    Do not download related files.
  --no-annotations      Do not download annotations.
  --no-auto-retry       Ask before retrying to download a file
  --retry-amount RETRY_AMOUNT
                        Number of times to retry a download
  --wait-time WAIT_TIME
                        Amount of seconds to wait before retrying
  --latest              Download latest version of a file if it exists
  --config FILE         Path to INI-type config file
  -m MANIFEST, --manifest MANIFEST
                        GDC download manifest file
```

- Download Manifest files.

```markdown

./gdc-client download -m gdc_manifest.2020-09-27-LUNG-miRNA-clinical.txt -d clinical/

./gdc-client download -m gdc_manifest.2020-09-27-LUNG-miRNA-isoform.txt -d isoform/

./gdc-client download -m gdc_manifest.2020-09-27-LUNG-miRNA-seq.txt -d miRNAseq/
```

- Check the files downloaded.

```markdown
$ cd clinical

# count total number of files downloaded
$ ls|wc
output
    522     522   19314

# count alive number of clinical information
$ grep -i vital_status */*xml|grep Alive|wc
output
    840     840   10920

# count dead number of clinical information
$ grep -i vital_status */*xml|grep -v Alive|wc
output
     298    3580  104444

# get sampes name/id:
$ grep -i vital_status */*xml|grep Alive |cut -d"." -f 3
output
TCGA-J2-8192
TCGA-J2-8192
TCGA-91-8499
TCGA-91-8499
TCGA-55-6986
TCGA-55-6986
TCGA-NJ-A4YG
TCGA-NJ-A4YG
......

# count total number of alive samples
$ grep -i vital_status */*xml|grep Alive |cut -d"." -f 3|sort -u|wc
output
    395     395    5135
```

## Data Wrangling

### Data of overall Information
  
```markdown
# count files number for clinical, isoform and miRNAseq
$ ls clinical|wc
output
     522     522   19314

$ ls isoform|wc
output
     567     567   20979

$ ls miRNAseq|wc
output
     567     567   20979
```

### Single sample

- Choose one file of clinical randomly to get the format of sample.

```markdown
$ ls clinical/57d733f6-2f8e-4f0c-b391-48d3cf7c1f87
output
nationwidechildrens.org_clinical.TCGA-75-7030.xml
```

- Follow the usage of [R - XML Files](https://www.tutorialspoint.com/r/r_xml_files.htm).

  - Open R

    The xml file is read by R using the function xmlParse(). It is stored as a list in R.

  - Reading XML File

    ```r
    # Load the package required to read XML files.
    library("XML")

    # Also load the other required package.
    library("methods")

    # Give the input file name to the function.
    result <- xmlParse(file = "nationwidechildrens.org_clinical.TCGA-75-7030.xml")

    # Print the result.
    print(result)
    ```

  - Get Number of Nodes Present in XML File

    ```r
    # Exract the root node form the xml file.
    rootnode <- xmlRoot(result)

    # Find number of nodes in the root.
    rootsize <- xmlSize(rootnode)

    # Print the result.
    print(rootsize)

    output
    [1] 2
    ```

  - Details of the First Node and second Node
  
    ```r
    # Exract the root node form the xml file.
    rootnode <- xmlRoot(result)

    # Print the result.
    print(rootnode[1])
    print(rootnode[2])
    ```

  - XML to Data Frame
  
    ```r
    # Convert the input xml file to a data frame.
    xmldataframe <- xmlToDataFrame(rootnode[2])
    print(xmldataframe)
    t(xmldataframe)
    write.table(t(xmldataframe),'tmp')
    ```

### Run R Scripts

- Complete scripts according to that of single sample.

- Run [Scripts](gdc.Rmd)
