
# coding: utf-8

# In[90]:

import itertools
import HTSeq
import difflib as df
import argparse
import os,sys,subprocess
import shutil,glob,re
import pprint
import ConfigParser,glob


# In[91]:

def setup_pipeline():
    config=ConfigParser.ConfigParser()
    path_to_config=os.path.join("pypeamplicon.cfg")
    config.read(path_to_config)
    sections=config.sections()
    r1=''
    r2=''
    se=''
    for section in sections:
        for software in config.options(section):
            if section == 'SOFTWARE':
                software_path=config.get(section,software)
                if os.path.isfile(software_path):
                    print "Path to {0} found".format(os.path.basename(software_path))
                    if "usearch" in os.path.basename(software_path):
                        usearch=software_path
                    elif "bbsplit" in os.path.basename(software_path):
                        bbsplit=software_path
                    elif "seqtk" in os.path.basename(software_path):
                        seqtk=software_path
                    elif "bbmap" in os.path.basename(software_path):
                        bbmap=software_path
                    elif "picard" in os.path.basename(software_path):
                        picard=software_path
                    elif "freebayes" in os.path.basename(software_path):
                        freebayes=software_path
                else:
                    raise Exception("Software {0} doesn't exist".format(os.path.basename(software_path)))
            elif section == "FQDIR":
                fqdir=config.get(section,software)
            elif section=="WORKINGDIR":
                cwd=config.get(section,software)
            elif section == "Ref":
                 ampref=config.get(section,software)
            elif section=="COMMANSUFFIX":
                suffix=config.get(section,software)
                if "_R1" in config.get(section,software):
                    r1=suffix
                elif "_R2"  in  config.get(section,software):
                    r2=suffix
                elif  "_SE" or " " in  config.get(section,software):
                    se=suffix
                else:
                    raise Exception("Please add correct suffix to fastq file")
            elif section == "CLUSTERFILTER":
                fcut=config.get(section,'filtercutoff')
                idcut=config.get(section,'idcutoff')
            else:
                raise Exception("Please correct the configeration file")


    return   usearch,bbsplit,seqtk,bbmap,picard,freebayes,fqdir,cwd,ampref,r1,r2,se,fcut,idcut  


# In[92]:

def merging_paired_reads(Read_Merged_OUT,R1_names,R2_names):
    print "Short read merging started....."
    merged_reads=[]
    try:
        print "Removing exsiting Read_Merged Directory"
        shutil.rmtree(Read_Merged_OUT)
        print "Creating Read_Merged Directory"
        os.mkdir(Read_Merged_OUT)
    except:
        print "Creating new Read_Merged Directory"
        os.mkdir(Read_Merged_OUT)

    for no,sample in enumerate(zip(R1_names,R2_names)):

        merged_file = os.path.join(Read_Merged_OUT,os.path.splitext(os.path.basename(sample[0]))[0] +"_merged.fastq")
        merged_reads.append(merged_file)
        Umerged_F_file = os.path.join(Read_Merged_OUT,os.path.splitext(os.path.basename(sample[0]))[0] +"_FUmerged.fastq")
        Umerged_R_file = os.path.join(Read_Merged_OUT,os.path.splitext(os.path.basename(sample[0]))[0] +"_RUmerged.fastq")
        report_file = os.path.join(Read_Merged_OUT,os.path.splitext(os.path.basename(sample[0]))[0] +"_report.txt")
        proc = subprocess.Popen([usearch, '-fastq_mergepairs', sample[0], '-reverse', sample[1], '-fastqout',merged_file, '-fastqout_notmerged_fwd',Umerged_F_file ,'-fastqout_notmerged_rev',Umerged_R_file,                             '-report',report_file, '-fastq_merge_maxee', '1','-fastq_minmergelen','100' ],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=False)
        #proc=subprocess.Popen(['seqtk'],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=False)
        (out, err) = proc.communicate()
        print "Program Output:",out
        print "Errors:", err
    if not err:
        print "Short read merging finished....."
        return merged_reads


# In[93]:

def read_files_from_location(qtrim_data, COMMANPREFIX_R1,COMMANPREFIX_R2,COMMANPREFIX_SE=''):
    if (COMMANPREFIX_R1 and COMMANPREFIX_R2):
        R1_file_names=sorted(glob.glob(os.path.join(qtrim_data,"*"+COMMANPREFIX_R1+"*")))
        R2_file_names=sorted(glob.glob(os.path.join(qtrim_data,"*"+COMMANPREFIX_R2+"*")))
    else:
        R1_file_names=[]
        R2_file_names=[]

    if COMMANPREFIX_SE:
        SE_file_names=glob.glob(os.path.join(qtrim_data,"*"+COMMANPREFIX_SE+"*"))
    else:
        SE_file_names=[]
    #files=os.listdir(qtrim_data)
    #for file in sorted(files):
    #	if COMMANPREFIX_R1 and COMMANPREFIX_R1 in file:
    #		print COMMANPREFIX_R1
    #		R1_file_names.append(os.path.join(qtrim_data,file))
    #	elif COMMANPREFIX_R2!='' and COMMANPREFIX_R2 in file:
    #		R2_file_names.append(os.path.join(qtrim_data,file))
    #	elif COMMANPREFIX_SE in file:
    #		SE_file_names.append(os.path.join(qtrim_data,file))
    return R1_file_names, R2_file_names,SE_file_names


# In[94]:

def read_clustering(Read_Merged_OUT,mode,destFolder,usearch,idcut):
    print "Short read clustering started....."
    read_clusters=destFolder
    merged_Fastq_dir=Read_Merged_OUT
    commanprefix_Merged="_merged.fastq"
    commanprefix_FU ="R1_FUmerged.fastq"
    commanprefix_RU ="R1_RUmerged.fastq"
    uR1_names,uR2_names,merged_reads = read_files_from_location(merged_Fastq_dir, commanprefix_FU,commanprefix_RU,commanprefix_Merged)
    for sample in merged_reads:
        cluster_file = os.path.join(read_clusters,os.path.splitext(os.path.basename(sample))[0] +"_cl")
        cns_file = os.path.join(read_clusters,os.path.splitext(os.path.basename(sample[0]))[0] +"_cns")
        proc = subprocess.Popen([usearch, '-cluster_fast', sample, '-id', idcut, '--clusters',cluster_file, '-consout',cns_file],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=False)
        (out, err) = proc.communicate()
        print "Program Output:",out
        print "Errors:", err
    print "Short read clustering finished...."



# In[95]:

def single_reads_prep(WORKING_DIR,SE_names):
    read_single=WORKING_DIR
    merged_reads=[]
    try:
        print "Removing exsiting Single-End Directory"
        shutil.rmtree(read_single)
        print "Creating New Single-End Directory"
        os.mkdir(read_single)
    except:
        print "Creating new Single-End Directory"
        os.mkdir(read_single)

    for read in SE_names:
        se_read=os.path.join(read_single,os.path.splitext(os.path.basename(read))[0] +"_merged.fastq")
        merged_reads.append(se_read)
        proc = subprocess.Popen(['ln', '-s', read, se_read],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=False)
        (out, err) = proc.communicate()
        print "Program Output:",out
        print "Errors:", err

    return merged_reads


# In[96]:

def filter_clusters(cluster_folder,filtered_clusters, fcut):
    print "Filtering clusters step has been started...."
    list_of_cluster=glob.glob(os.path.join(cluster_folder,"*_cl*"))
    no_intial_clusters=len(list_of_cluster)
    for cluster in list_of_cluster:
        proc = subprocess.Popen(['grep', '-c', ">", cluster ],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=False)
        (out, err) = proc.communicate()
        print "Program Output:",out
        if not err:
            if int(out) > int(fcut):
                shutil.move(cluster, os.path.join(filtered_clusters,os.path.basename(cluster)))

    no_filterd_clusters=len(glob.glob(os.path.join(filtered_clusters,"*_cl*")))
    print "Initial Clusters : {0}".format(no_intial_clusters)
    print "No of clusters filtered :{0}".format(no_filterd_clusters)
    print "Filtering clusters step has been finished...."


# In[97]:

def spliting_referance(referance, destination):
    reads = HTSeq.FastaReader( AMPLICON_FASTA )


    list_of_amplicon=[]
    for read in reads:
        fasta_file_name = os.path.join(destination,read.name.strip()+".fa")
        list_of_amplicon.append(fasta_file_name)
        with open(fasta_file_name,"w") as f:
            read.write_to_fasta_file(f)


    bb_list_of_amplicon=",".join(list_of_amplicon)
    return bb_list_of_amplicon


# In[98]:

def making_directory(path, dirname):
    path_to_dir=path
    try:
        print "Removing exsiting {0} Directory".format(dirname)
        shutil.rmtree(path_to_dir)
        print "Creating New {0} Directory".format(dirname)
        os.mkdir(path_to_dir)
    except:
        print "Creating new {0} Directory".format(dirname)
        os.mkdir(path_to_dir)


# In[99]:

def populating_directories(WORKING_DIR):
    cluster_folder=os.path.join(WORKING_DIR,"Read-Clusters")
    align_clusters=os.path.join(WORKING_DIR,"Align-Clusters")
    filtered_clusters=os.path.join(cluster_folder,"Filtered-Clusters")
    making_directory(cluster_folder, "Read-Clusters")
    making_directory(align_clusters, "Align-Clusters")
    making_directory(filtered_clusters, "Filtered-Clusters")

    return cluster_folder,align_clusters,filtered_clusters


# In[100]:

def spliting_fastq_per_referance(cluster_folder,filtered_clusters,align_clusters,bb_list_of_amplicon,bbsplit):
    print "Spliting fastq to clusters step has been started...."

    filtered_clusters=os.path.join(cluster_folder,"Filtered-Clusters")
    #bb_list_of_amplicon=",".join(list_of_amplicon)

    for cluster in glob.glob(os.path.join(filtered_clusters,"*_cl*")):
        algined_cluster_dir=os.path.join(align_clusters,os.path.basename(cluster))
        try:
            shutil.rmtree(algined_cluster_dir)
            os.mkdir(algined_cluster_dir)
        except:
            os.mkdir(algined_cluster_dir)
        proc = subprocess.Popen([bbsplit,  'ref='+bb_list_of_amplicon, 'in='+cluster,re.escape('basename=')+r'o%'+re.escape('.fq')],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=False)
        (out, err) = proc.communicate()
        print "Program Output:",out
        print "Errors:", err
        for f in glob.glob('o*.fq'):
            if os.path.getsize(f)>0:
                shutil.move(f, algined_cluster_dir)
            else:
                 os.remove(f)
    print "Spliting fastq to clusters step has been finished...."


# In[101]:

def converting_fasta_to_fastq(align_clusters,merged_reads,seqtk):
    print "Converting Fasta to Fastq step has been started...."

    if os.path.isfile(os.path.join(align_clusters,"temp.fastq")):
        shutil.rmtree(os.path.join(align_clusters,"temp.fastq"))
        merged_fq=open(os.path.join(align_clusters,"temp.fastq"), "a")
    else:
        merged_fq=open(os.path.join(align_clusters,"temp.fastq"), "a")

    print merged_reads
    for fq  in merged_reads:
        fq_reads = HTSeq.FastqReader( fq )
        for r in fq_reads:
            r.write_to_fastq_file(merged_fq)
    merged_fq.close()

    fastq_dir=glob.glob(os.path.join(align_clusters,'*_cl*','*.fq'))

    reads_merged = os.path.join(align_clusters,"temp.fastq")
    for fq in fastq_dir:
        cl_files= HTSeq.FastqReader(fq)
        reads_merged_header_l = os.path.join(align_clusters,"temp.list")
        proc1 = subprocess.Popen(['grep', '^@M', fq],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=False)
        (out1, err1) = proc1.communicate()
        if not err1:
                with open(reads_merged_header_l,"w") as f:
                    f.write(out1.replace("@",""))

        else:
            print "Errors:", err1

        proc2 = subprocess.Popen([seqtk, 'subseq', reads_merged, reads_merged_header_l],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=False)
        (out2, err2) = proc2.communicate()
        if not err2:
            new_fq_w_q=os.path.splitext(fq)[0]+"_q.fq"
            with open( new_fq_w_q,"w") as f:
                f.write(out2)


    print "Converting Fasta to Fastq step has been finished...."


# In[142]:

def mapping_to_reference(align_clusters,AMPLICON_FASTA,bbmap,picard):
    print "Mapping  Fastq to amplicons  step has been started...."
    fastq_dir=glob.glob(os.path.join(align_clusters,'*_cl*','*_q.fq'))
    for fq in fastq_dir:
        name=os.path.splitext(os.path.basename(fq))[0]

        sam_file_dir=os.path.join(os.path.splitext(fq)[0]+".sam")
        bam_file_dir_sorted=os.path.join(os.path.splitext(fq)[0]+"_sorted.bam")
        proc = subprocess.Popen([bbmap,  'ref='+AMPLICON_FASTA, 'in='+fq,'out='+sam_file_dir],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=False)
        (out, err) = proc.communicate()
        proc2 = subprocess.Popen(['java', '-jar' ,'-Xmx4g', picard, 'SamFormatConverter','I='+sam_file_dir, 'OUTPUT='+os.path.join('/dev/stdout')],stdout=subprocess.PIPE,shell=False)
        #(out2, err2) = proc2.communicate()
        #print "Program Output:",out2
        #print "Errors:", err2
        sid=['ID={0}'.format("".join(os.path.splitext(fq)[0].split("/")[-2].split('_merged_')[0].split("_")[0]))]
        lb=['LB={0}'.format( os.path.splitext(fq)[0].split("/")[-2].split('_merged_')[0])]
        pl=['Pl=illumina']
        pu=['PU=None']
        sm=['SM={0}'.format("".join(os.path.splitext(fq)[0].split("/")[-2].split('_merged_')[0].split("_")[0]))]
        proc3 = subprocess.Popen(['java', '-jar' ,'-Xmx4g', picard,'SortSam','VALIDATION_STRINGENCY=LENIENT', 'INPUT='+os.path.join('/dev/stdin'), 'OUTPUT='+os.path.join('/dev/stdout'),'SORT_ORDER=coordinate'],stdin=proc2.stdout,stdout=subprocess.PIPE,shell=False)
        #(out1, err1) = proc1.communicate()
        #print "Errors:", err1

        proc4 = subprocess.Popen(['java', '-jar' ,'-Xmx4g', picard, 'AddOrReplaceReadGroups','I='+os.path.join('/dev/stdin'), 'OUTPUT='+bam_file_dir_sorted]+sid+lb+pl+pu+sm,stdout=subprocess.PIPE, stderr=subprocess.STDOUT,stdin=proc3.stdout,shell=False)
        (out4, err4) = proc4.communicate()
        print "Program Output:",out4
        print "Errors:", err4


    print "Mapping  Fastq to amplicons  step has been finished...."


# In[130]:

def moving_sam_to_bam(align_clusters, WORKING_DIR,picard):
    print "Moving Sorted Bam files...."
    bam_list=[]


    bam_f_sorted=glob.glob(os.path.join(align_clusters,'*_cl*','*_q_sorted.bam'))
    print bam_f_sorted
    amplicon_dict=dict() 
    for bam in bam_f_sorted:
        name=str(os.path.splitext(os.path.basename(bam))[0])
        cluster=os.path.splitext(bam)[0].split("/")[-2].replace('_R1_merged',"")
        key=name.lstrip('o').rstrip('_q_sorted')
        if key in amplicon_dict.keys():
            new_folder=os.path.join(align_clusters,key)
            shutil.copy(bam,os.path.join(new_folder,key+cluster+"_sorted.bam"))
            amplicon_dict[key].append(os.path.join(new_folder,key+cluster+"_sorted.bam"))
        else:
            new_folder=os.path.join(align_clusters,key)
            try:
                os.mkdir(new_folder)
            except:
                shutil.rmtree(new_folder)
                os.mkdir(new_folder)

            shutil.copy(bam,os.path.join(new_folder,key+cluster+"_sorted.bam"))
            amplicon_dict[key]=[os.path.join(new_folder,key+cluster+"_sorted.bam")]



    print "Adding Groups to Sam Files and Coverting Sam to Bam...."

    for key,path in amplicon_dict.items():
        sorted_bam_list=["INPUT="+os.path.join(x.replace("./","")) for x in path]
        out_dir=os.path.dirname(path[0])
        name_merged_bam=os.path.splitext(path[0])[0].split('/')[-2]
        out_dir_merged_bam=os.path.join(out_dir,name_merged_bam+"_merged.bam")
        proc = subprocess.Popen(['java', '-jar' ,'-Xmx4g', picard,'MergeSamFiles','USE_THREADING=true','VALIDATION_STRINGENCY=LENIENT','OUTPUT='+out_dir_merged_bam]+sorted_bam_list,stdout=subprocess.PIPE,shell=False)
        (out, err) = proc.communicate()
        print "Program Output:",out
        print "Errors:", err
        #proc2 = subprocess.Popen(['java', '-jar' ,'-Xmx4g', picard, 'SamFormatConverter','I='+os.path.join('/dev/stdin'), 'OUTPUT='+out_dir_merged_bam],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,stdin=proc.stdout,shell=False)
        #(out2, err2) = proc2.communicate()
        #print "Program Output:",out2
        #print "Errors:", err2
        proc3 = subprocess.Popen(['java', '-jar' ,'-Xmx4g', picard, 'BuildBamIndex','I='+out_dir_merged_bam, 'O='+out_dir_merged_bam+".bai"],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=False)
        (out3, err3) = proc3.communicate()
        print "Program Output:",out3
        print "Errors:", err3
        bam_list.append(out_dir_merged_bam)

    print "Creating final Merged Bamfile"

    sorted_bam_list=["INPUT="+os.path.join(x.replace("./","")) for x in bam_list]
    proc4 = subprocess.Popen(['java', '-jar' ,'-Xmx4g', picard,'MergeSamFiles','USE_THREADING=true','VALIDATION_STRINGENCY=LENIENT','OUTPUT='+os.path.join(WORKING_DIR, "Final_Merged.bam")]+sorted_bam_list,stdout=subprocess.PIPE,shell=False)
    (out4, err4) = proc4.communicate()
    print "Program Output:",out4
    print "Errors:", err4

    proc6 = subprocess.Popen(['java', '-jar' ,'-Xmx4g', picard,'SortSam','VALIDATION_STRINGENCY=LENIENT', 'INPUT='+os.path.join(WORKING_DIR, "Final_Merged.bam"), 'OUTPUT='+os.path.join(WORKING_DIR, "Final_Merged_csorted.bam"),'SORT_ORDER=coordinate'],stdout=subprocess.PIPE,shell=False)
    (out6, err6) = proc6.communicate()
    print "Errors:", err6


    proc5 = subprocess.Popen(['java', '-jar' ,'-Xmx4g', picard, 'BuildBamIndex','I='+os.path.join(WORKING_DIR, "Final_Merged_csorted.bam"), 'O='+os.path.join(WORKING_DIR, "Final_Merged_csorted.bam")+".bai"],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=False)
    (out5, err5) = proc5.communicate()
    print "Program Output:",out5
    print "Errors:", err5


# In[132]:

def generate_vcf(WORKING_DIR, AMPLICON_FASTA, freebayes):
    proc = subprocess.Popen([freebayes, '-f' , AMPLICON_FASTA, os.path.join(WORKING_DIR, "Final_Merged_csorted.bam"), '-v',os.path.join(WORKING_DIR, "Final_Merged.vcf")],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=False)
    (out, err) = proc.communicate()
    print "Program Output:",out
    print "Errors:", err

    print "Analysis is done...."


if __name__ == '__main__':
    mode=''
    usearch,bbsplit,seqtk,bbmap,picard,freebayes,fqdir,cwd,ampref,r1,r2,se,fcut,idcut=setup_pipeline()
    FASTQ_DIR=os.path.join(fqdir)
    WORKING_DIR=os.path.join(cwd)

    single_read_loc=''
    AMPLICON_FASTA = os.path.join(WORKING_DIR,ampref)
    COMMANPREFIX_R1=r1 #Check your Q-Trrimmed filenames
    COMMANPREFIX_R2=r2 #Check your Q-Trrimmed filenames
    COMMANPREFIX_SE=se #Check your Q-Trrimmed filenames
    cluster_folder,align_clusters,filtered_clusters=populating_directories(WORKING_DIR)


    R1_names,R2_names,SE_names = read_files_from_location(FASTQ_DIR, COMMANPREFIX_R1,COMMANPREFIX_R2,COMMANPREFIX_SE)
    if R1_names and R2_names:
        print "Paired-end data has been detected. Switching to paired end mode......"
        mode='p'
    else:
        print "Single-end data has been detected. Skipping read merging step......"
        mode='s'

    if mode=='p':
        loc=os.path.join(WORKING_DIR,"Read-Merged")
        merged_reads=merging_paired_reads(loc,R1_names,R2_names)
        print merged_reads

    if mode=='s':
        loc=os.path.join(WORKING_DIR,"Single-End")
        merged_reads=single_reads_prep(loc,SE_names)


    # In[136]:


    read_clustering(loc,mode,cluster_folder,usearch,idcut)


    # In[137]:

    filter_clusters(cluster_folder,filtered_clusters, fcut)


    # In[138]:

    bb_list_of_amplicon=spliting_referance(AMPLICON_FASTA, align_clusters)


    # In[139]:

    spliting_fastq_per_referance(cluster_folder,filtered_clusters,align_clusters,bb_list_of_amplicon,bbsplit)


    # In[140]:

    converting_fasta_to_fastq(align_clusters,merged_reads, seqtk)


    # In[143]:

    mapping_to_reference(align_clusters, AMPLICON_FASTA,bbmap,picard)


    # In[144]:

    moving_sam_to_bam(align_clusters, WORKING_DIR,picard)


    # In[145]:

    generate_vcf(WORKING_DIR, AMPLICON_FASTA,freebayes)

