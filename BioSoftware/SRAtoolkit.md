# SRAtoolkit Install

1. [Download](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) corresponding version. Upload the install package to the server. 

2. Unzip: `tar zvxf ~/bioinformatics/software/sratoolkit.3.0.0-ubuntu64.tar.gz -C ~/bioinformatics/software`

3. Configuration: `~/bioinformatics/software/sratoolkit.3.0.0-ubuntu64/bin/vdb-config --interactive`
 
   - Create a new folder called "sratoolkit" to save user information 

   - After entering the terminal，type "C" to choose CACHE. The location of user-repository is the location of new folder "sratoolkit". The process-local location is the location of sratoolkit.3.0.0-ubuntu64 folder. 
  
   - Type "A" to choose AWS. Choose "report cloud instance identity".
   
   - Type "S" to save. Type "X" to exit.
  
4. Test installation: `~/bioinformatics/software/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump -h`

5. Set environmental variables：  
    ```
    echo 'export PATH=~/bioinformatics/software/sratoolkit.3.0.0-ubuntu64/bin:$PATH'  >> ~/.bashrc 
    soruce ~/.bashrc
    ```
6. Check environmental variables setting: `fastq-dump`

# Download SRA data in batch

1. Search dataset by GEO number, such as [GSE153481](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153481).

2. Click "SRA Run Selector".

3. Click  "Accession List" to download file list. Upload the file list to the server. 

4. Download in batch: `prefetch --option-file SRR_Acc_List.txt`
