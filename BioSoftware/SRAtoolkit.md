# SRAtoolkit Install

1. [Download](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) corresponding version. Upload the install package to the server. 

2. Unzip: `tar zvxf ~/bioinformatics/software/sratoolkit.3.0.0-ubuntu64.tar.gz -C ~/bioinformatics/software`

3. Configuration: `~/bioinformatics/software/sratoolkit.3.0.0-ubuntu64/bin/vdb-config --interactive`
 
   - Create a new folder called "sratoolkit" to save user information 

   - After entering the terminal，type "C" to choose CACHE. The location of user-repository is the location of new folder "sratoolkit". The process-local location is the location of sratoolkit.3.0.0-ubuntu64 folder. 
  
   - Type "A" to choose AWS. Choose "report cloud instance identity".
   
   - Type "S" to save. Type "X" to exit.
  
4. Test installation: `~/bioinformatics/software/sratoolkit.3.0.0-ubuntu64/bin/vdb-config --interactive`

7. 设置环境变量：  
    ```
    echo 'export PATH=~/bioinformatics/software/sratoolkit.3.0.0-ubuntu64/bin:$PATH'  >> ~/.bashrc 
    soruce ~/.bashrc
    ```
6. 检查环境变量设置是否成功: `fastq-dump`

# 批量下载SRA文件

1. 根据GEO号查找相关数据集，例如 [GSE153481](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153481 )
2. 点击 "SRA Run Selector"
3. 点击 "Accession List" 下载文件列表
4. 将下载的文件通过WinSCP上传到服务器相应位置
5. 批量下载 `prefetch --option-file SRR_Acc_List.txt`
