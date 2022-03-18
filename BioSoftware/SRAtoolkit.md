# 安装SRAtoolkit

1. 从[官网](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)下载压缩包，并上传到服务器
2. 解压压缩包 `tar zvxf ~/bioinformatics/software/sratoolkit.3.0.0-ubuntu64.tar.gz -C ~/bioinformatics/software`
3. 运行配置环境`~/bioinformatics/software/sratoolkit.3.0.0-ubuntu64/bin/vdb-config --interactive`
   - 新建空文件夹 sratoolkit 用于保存用户信息
   - 进入终端后，按C选择CACHE，location of user-repository 选择新建的空文件夹sratoolkit，process-local location选择安装包解压成的文件夹sratoolkit.3.0.0-ubuntu64
   - 按A选择AWS，选择report cloud instance identity
   - 按S保存设置，按X退出
4. 执行 `~/bioinformatics/software/sratoolkit.3.0.0-ubuntu64/bin/vdb-config --interactive`   检查安装结果

