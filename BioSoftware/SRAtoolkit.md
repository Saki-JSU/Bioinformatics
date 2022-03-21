# 安装SRAtoolkit

1. 从[官网](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)下载压缩包，并上传到服务器
2. 解压压缩包: `tar zvxf ~/bioinformatics/software/sratoolkit.3.0.0-ubuntu64.tar.gz -C ~/bioinformatics/software`
3. 运行配置环境: `~/bioinformatics/software/sratoolkit.3.0.0-ubuntu64/bin/vdb-config --interactive`
   - 新建空文件夹 sratoolkit 用于保存用户信息
   - 进入终端后，按C选择CACHE，location of user-repository 选择新建的空文件夹sratoolkit，process-local location选择安装包解压成的文件夹sratoolkit.3.0.0-ubuntu64
   - 按A选择AWS，选择report cloud instance identity
   - 按S保存设置，按X退出
4. 执行 `~/bioinformatics/software/sratoolkit.3.0.0-ubuntu64/bin/vdb-config --interactive`   检查安装结果
5. 设置环境变量：  
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
