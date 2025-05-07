$HOSTNAME = ""
params.outdir = 'results'  

def pathChecker(input, path, type){
	def cmd = "mkdir -p check && mv ${input} check/. "
	if (!input || input.empty()){
		input = file(path).getName().toString()
		cmd = "mkdir -p check && cd check && ln -s ${path} ${input} && cd .."
		if (path.indexOf('s3:') > -1 || path.indexOf('S3:') >-1){
			def recursive = (type == "folder") ? "--recursive" : ""
			cmd = "mkdir -p check && cd check && aws s3 cp ${recursive} ${path} ${workDir}/${input} && ln -s ${workDir}/${input} . && cd .."
		} else if (path.indexOf('gs:') > -1 || path.indexOf('GS:') >-1){
			if (type == "folder"){
				cmd = "mkdir -p check ${workDir}/${input} && cd check && gsutil rsync -r ${path} ${workDir}/${input} && cp -R ${workDir}/${input} . && cd .."
			} else {
				cmd = "mkdir -p check && cd check && gsutil cp ${path} ${workDir}/${input} && cp -R ${workDir}/${input} . && cd .."
			}
		} else if (path.indexOf('/') == -1){
			cmd = ""
		}
}
	return [cmd,input]
}
if (!params.groups_file){params.groups_file = ""} 
if (!params.compare_file){params.compare_file = ""} 
if (!params.gtf){params.gtf = ""} 
if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.run_DESeq2){params.run_DESeq2 = ""} 
if (!params.genome){params.genome = ""} 
if (!params.run_FeatureCounts){params.run_FeatureCounts = ""} 
if (!params.run_limmaVoom){params.run_limmaVoom = ""} 
if (!params.run_gsea_DESeq2){params.run_gsea_DESeq2 = ""} 
if (!params.run_gsea_LimmaVoom){params.run_gsea_LimmaVoom = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)

g_3_1_g37_24 = params.groups_file && file(params.groups_file, type: 'any').exists() ? file(params.groups_file, type: 'any') : ch_empty_file_1
g_3_1_g37_42 = params.groups_file && file(params.groups_file, type: 'any').exists() ? file(params.groups_file, type: 'any') : ch_empty_file_1
g_4_2_g37_24 = params.compare_file && file(params.compare_file, type: 'any').exists() ? file(params.compare_file, type: 'any') : ch_empty_file_2
g_4_2_g37_42 = params.compare_file && file(params.compare_file, type: 'any').exists() ? file(params.compare_file, type: 'any') : ch_empty_file_2
g_12_3_g_27 = file(params.gtf, type: 'any')
g_12_1_g14_21 = file(params.gtf, type: 'any')
if (params.reads){
Channel
	.fromFilePairs( params.reads,checkExists:true , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 ) 
	.set{g_15_0_g_31}
 (g_15_0_g_33) = [g_15_0_g_31]
 } else {  
	g_15_0_g_31 = Channel.empty()
	g_15_0_g_33 = Channel.empty()
 }

Channel.value(params.mate).set{g_16_1_g_31}
(g_16_1_g_33,g_16_1_g_27,g_16_1_g14_30,g_16_0_g14_31) = [g_16_1_g_31,g_16_1_g_31,g_16_1_g_31,g_16_1_g_31]
Channel.value(params.run_DESeq2).set{g_21_3_g37_24}
g_30_0_g14_21 = file(params.genome, type: 'any')
Channel.value(params.run_FeatureCounts).set{g_36_0_g_35}
Channel.value(params.run_limmaVoom).set{g_38_3_g37_42}
Channel.value(params.run_gsea_DESeq2).set{g_39_2_g37_33}
Channel.value(params.run_gsea_LimmaVoom).set{g_40_2_g37_41}

//* params.run_Adapter_Removal =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Adapter_Removal"
//* @style @multicolumn:{seed_mismatches, palindrome_clip_threshold, simple_clip_threshold} @condition:{Tool_for_Adapter_Removal="trimmomatic", seed_mismatches, palindrome_clip_threshold, simple_clip_threshold}, {Tool_for_Adapter_Removal="fastx_clipper", discard_non_clipped}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 5
    $MEMORY = 8
}
//* platform
//* platform
//* autofill

process Adapter_Removal {

input:
 tuple val(name), file(reads)
 val mate

output:
 tuple val(name), file("reads/*.fastq.gz")  ,emit:g_31_reads01_g14_31 
 path "*.{fastx,trimmomatic}.log"  ,emit:g_31_log_file11 

container 'quay.io/viascientific/trimmomatic:1.0'

when:
(params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal

shell:
phred = params.Adapter_Removal.phred
Tool_for_Adapter_Removal = params.Adapter_Removal.Tool_for_Adapter_Removal
Adapter_Sequence = params.Adapter_Removal.Adapter_Sequence
//trimmomatic_inputs
min_length = params.Adapter_Removal.min_length
seed_mismatches = params.Adapter_Removal.seed_mismatches
palindrome_clip_threshold = params.Adapter_Removal.palindrome_clip_threshold
simple_clip_threshold = params.Adapter_Removal.simple_clip_threshold

//fastx_clipper_inputs
discard_non_clipped = params.Adapter_Removal.discard_non_clipped
    
remove_previous_reads = params.Adapter_Removal.remove_previous_reads
discard_non_clipped_text = ""
if (discard_non_clipped == "yes") {discard_non_clipped_text = "-c"}
nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] 
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1] }
} 
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage;
 use Cwd qw();
 
runCmd("mkdir reads adapter unpaired");

open(OUT, ">adapter/adapter.fa");
my @adaps=split(/\n/,"!{Adapter_Sequence}");
my $i=1;
foreach my $adap (@adaps)
{
 print OUT ">adapter$i\\n$adap\\n";
 $i++;
}
close(OUT);

my $quality="!{phred}";
print "fastq quality: $quality\\n";
print "tool: !{Tool_for_Adapter_Removal}\\n";

if ("!{mate}" eq "pair") {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic PE -threads !{task.cpus} -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq.gz unpaired/!{name}.1.fastq.unpaired.gz reads/!{name}.2.fastq.gz unpaired/!{name}.2.fastq.unpaired.gz ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        print "Fastx_clipper is not suitable for paired reads.";
    }
} else {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic SE -threads !{task.cpus}  -phred${quality} !{file1} reads/!{name}.fastq.gz ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        runCmd("fastx_clipper  -Q $quality -a !{Adapter_Sequence} -l !{min_length} !{discard_non_clipped_text} -v -i !{file1} -o reads/!{name}.fastq.gz > !{name}.fastx.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir = join '/', @paths;
    $inputsdir .= "/work";
    print "INFO: inputs reads will be removed if they are located in the $workdir $inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            runCmd("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}


##Subroutines
sub runCmd {
    my ($com) = @_;
    if ($com eq ""){
		return "";
    }
    my $error = system(@_);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}

//* params.run_FastQC =  "no"  //* @dropdown @options:"yes","no"
if (params.run_FastQC == "no") { println "INFO: FastQC will be skipped"}


process FastQC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.(html|zip)$/) "fastqc_report/$filename"}
input:
 tuple val(name), file(reads)
 val mate

output:
 path '*.{html,zip}'  ,emit:g_33_FastQCout00 

errorStrategy 'retry'
maxRetries 5

script:
nameAll = reads.toString()
if (nameAll.contains('.gz')) {
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    file =  nameAll 
    runGzip = ''
}
"""
if [ "${params.run_FastQC}" == "yes" ]; then
    ${runGzip}
    fastqc ${file} 
else
    touch process.skiped.html
fi
"""
}

//* @style @array:{run_name,run_parameters} @multicolumn:{run_name,run_parameters}

process featureCounts_Prep {

input:
 val run_featureCounts

output:
 val run_params  ,emit:g_35_run_parameters02_g_27 

when:
run_featureCounts == "yes"

script:
run_name = params.featureCounts_Prep.run_name
run_parameters = params.featureCounts_Prep.run_parameters
sense_antisense = params.featureCounts_Prep.sense_antisense

//define run_name and run_parameters in map item and push into run_params array
run_params = []
for (i = 0; i < run_parameters.size(); i++) {
   map = [:]
   map["run_name"] = run_name[i].replaceAll(" ","_").replaceAll(",","_").replaceAll(";","_").replaceAll("'","_").replaceAll('"',"_")
   map["run_parameters"] = run_parameters[i]
   run_params[i] = map
}
templateRunParams = run_parameters[0] ? run_parameters[0] : "-g gene_id -s 0 -Q 20 -T 2 -B -d 50 -D 1000 -C --fracOverlap 0 --minOverlap 1"
if (sense_antisense == "Yes"){
   map = [:]
   map["run_name"] = "gene_id_forward_expression"
   map["run_parameters"] = templateRunParams.replaceAll("transcript_id","gene_id").replaceAll("-s 0","-s 1").replaceAll("-s 2","-s 1")
   run_params.push(map)
   map = [:]
   map["run_name"] = "gene_id_reverse_expression"
   map["run_parameters"] = templateRunParams.replaceAll("transcript_id","gene_id").replaceAll("-s 0","-s 2").replaceAll("-s 1","-s 2")
   run_params.push(map)
   map = [:]
   map["run_name"] = "transcript_id_forward_expression"
   map["run_parameters"] = templateRunParams.replaceAll("gene_id","transcript_id").replaceAll("-s 0","-s 1").replaceAll("-s 2","-s 1")
   run_params.push(map)
   map = [:]
   map["run_name"] = "transcript_id_reverse_expression"
   map["run_parameters"] = templateRunParams.replaceAll("gene_id","transcript_id").replaceAll("-s 0","-s 2").replaceAll("-s 1","-s 2")
   run_params.push(map)
}
"""
"""

}

build_STAR_index = params.STAR_Module_Check_Build_STAR_Index.build_STAR_index

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 5
    $MEMORY = 150
}
//* platform
//* platform
//* autofill

process STAR_Module_Check_Build_STAR_Index {

input:
 path genome
 path gtf

output:
 path "STARIndex"  ,emit:g14_21_starIndex00_g14_26 

when:
build_STAR_index == true && ((params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR)

script:
star_build_parameters = params.STAR_Module_Check_Build_STAR_Index.star_build_parameters
newDirName = "STARIndex" 
"""
if [ ! -e "${params.star_index}/SA" ] ; then
    echo "STAR index not found"
    mkdir -p $newDirName 
    STAR --runMode genomeGenerate ${star_build_parameters} --genomeDir $newDirName --genomeFastaFiles ${genome} --sjdbGTFfile ${gtf}
else 
	ln -s ${params.star_index} STARIndex
fi

"""





}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 10
}
//* platform
//* platform
//* autofill

process STAR_Module_check_STAR_files {

input:
 path star

output:
 path "*/${star2}" ,optional:true  ,emit:g14_26_starIndex02_g14_31 

container 'quay.io/viascientific/pipeline_base_image:1.0'
stageInMode 'copy'

when:
(params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR

script:
(cmd, star2) = pathChecker(star, params.star_index, "folder")
"""
$cmd
"""
}

//* params.star_index =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 10
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process STAR_Module_Map_STAR {

input:
 val mate
 tuple val(name), file(reads)
 path star_index

output:
 tuple val(name), file("${name}Log.final.out")  ,emit:g14_31_outputFileOut00_g14_18 
 tuple val(name), file("${name}.flagstat.txt")  ,emit:g14_31_outputFileTxt11 
 tuple val(name), file("${name}Log.out")  ,emit:g14_31_logOut22 
 tuple val(name), file("${name}.bam")  ,emit:g14_31_mapped_reads30_g14_30 
 tuple val(name), file("${name}SJ.out.tab")  ,emit:g14_31_outputFileTab44 
 tuple val(name), file("${name}Log.progress.out")  ,emit:g14_31_progressOut55 
 tuple val(name), file("${name}Aligned.toTranscriptome.out.bam") ,optional:true  ,emit:g14_31_transcriptome_bam60_g14_15 

when:
(params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR

script:
params_STAR = params.STAR_Module_Map_STAR.params_STAR
sense_antisense = params.STAR_Module_Map_STAR.sense_antisense
transcriptomeSAM = ""
if (params.run_Salmon_after_STAR && params.run_Salmon_after_STAR == "yes" && params_STAR.indexOf("--quantMode") < 0){
	transcriptomeSAM = " --quantMode TranscriptomeSAM "
}

"""
STAR ${params_STAR} ${transcriptomeSAM} --genomeDir ${star_index} --readFilesCommand zcat --readFilesIn $reads --outFileNamePrefix ${name}
echo "Alignment completed."
if [ ! -e "${name}Aligned.toTranscriptome.out.bam" -a -e "${name}Aligned.toTranscriptome.out.sam" ] ; then
    samtools view -S -b ${name}Aligned.toTranscriptome.out.sam > ${name}Aligned.toTranscriptome.out.bam
elif [ ! -e "${name}Aligned.out.bam" -a -e "${name}Aligned.out.sam" ] ; then
    samtools view -S -b ${name}Aligned.out.sam > ${name}Aligned.out.bam
fi
rm -rf *.sam
if [ -e "${name}Aligned.sortedByCoord.out.bam" ] ; then
    mv ${name}Aligned.sortedByCoord.out.bam ${name}.bam
elif [ -e "${name}Aligned.out.bam" ] ; then
    mv ${name}Aligned.out.bam ${name}.bam
fi

samtools flagstat ${name}.bam > ${name}.flagstat.txt
"""


}


process STAR_Module_STAR_Summary {

input:
 tuple val(name), file(alignSum)

output:
 path "*.tsv"  ,emit:g14_18_outputFile00_g14_11 
 val "star_alignment_sum"  ,emit:g14_18_name11_g14_11 

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my %tsv;
my @headers = ();
my $name = '!{name}';


alteredAligned();

my @keys = keys %tsv;
my $summary = "$name"."_star_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
	my $values = join("\\t", @{ $tsv{$key} });
	`echo "$values" >> $summary`;
}


sub alteredAligned
{
	my @files = qw(!{alignSum});
	my $multimappedSum;
	my $alignedSum;
	my $inputCountSum;
	push(@headers, "Sample");
    push(@headers, "Total Reads");
	push(@headers, "Multimapped Reads Aligned (STAR)");
	push(@headers, "Unique Reads Aligned (STAR)");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($inputCount = `cat $file | grep 'Number of input reads' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($aligned = `cat $file | grep 'Uniquely mapped reads number' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($multimapped = `cat $file | grep 'Number of reads mapped to multiple loci' | awk '{sum+=\\$9} END {print sum}'`);
		$multimappedSum += int($multimapped);
        $alignedSum += int($aligned);
        $inputCountSum += int($inputCount);
	}
	$tsv{$name} = [$name, $inputCountSum];
	push(@{$tsv{$name}}, $multimappedSum);
	push(@{$tsv{$name}}, $alignedSum);
}

sub runCommand {
    my ($com) = @_;
    if ($com eq ""){
		return "";
    }
    my $error = system(@_);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}


process STAR_Module_merge_tsv_files_with_same_header {

input:
 path tsv
 val outputFileName

output:
 path "${name}.tsv"  ,emit:g14_11_outputFileTSV00 

errorStrategy 'retry'
maxRetries 3

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 10
}
//* platform
//* platform
//* autofill

process STAR_Module_Merge_Bam_and_create_sense_antisense {

input:
 tuple val(oldname), file(bamfiles)
 val mate

output:
 path "*_sorted.bam.bai"  ,emit:g14_30_bam_index00 
 tuple val(oldname), file("*_sorted.bam")  ,emit:g14_30_bamFile10_g_27 

errorStrategy 'retry'
maxRetries 2

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
    samtools merge !{oldname}.bam !{bamfiles.join(" ")} && samtools sort -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
else
    mv !{bamfiles.join(" ")} !{oldname}.bam 2>/dev/null || true
    samtools sort  -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
fi


'''
}

//* params.gtf =  ""  //* @input


process featureCounts {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*$/) "Counts/$filename"}
input:
 tuple val(name), file(bam)
 val paired
 each run_params
 path gtf

output:
 path "*"  ,emit:g_27_outputFileTSV00_g_28 

container 'quay.io/biocontainers/subread:1.6.4--h84994c4_1'

script:
pairText = ""
if (paired == "pair"){
    pairText = "-p"
}

run_name = run_params["run_name"] 
run_parameters = run_params["run_parameters"] 

"""
featureCounts ${pairText} ${run_parameters} -a ${gtf} -o ${name}@${run_name}@fCounts.txt ${bam}
## remove first line
sed -i '1d' ${name}@${run_name}@fCounts.txt

"""
}

//* autofill
//* platform
//* platform
//* autofill

process featureCounts_summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_featureCounts.tsv$/) "Counts_Summary/$filename"}
input:
 path featureCountsOut

output:
 path "*_featureCounts.tsv"  ,emit:g_28_outputFile00_g37_24 
 path "*_featureCounts.sum.tsv"  ,emit:g_28_outFileTSV11 

shell:
'''
#!/usr/bin/env perl

# Step 1: Merge count files
my %tf = ( expected_count => 6 );
my @run_name=();
chomp(my $contents = `ls *@fCounts.txt`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
    $file=~/(.*)\\@(.*)\\@fCounts\\.txt/;
    my $runname = $2;
    push(@run_name, $runname) unless grep{$_ eq $runname} @run_name;
}


my @expectedCount_ar = ("expected_count");
for($l = 0; $l <= $#run_name; $l++) {
    my $runName = $run_name[$l];
    for($ll = 0; $ll <= $#expectedCount_ar; $ll++) {
        my $expectedCount = $expectedCount_ar[$ll];
    
        my @a=();
        my %b=();
        my %c=();
        my $i=0;
        chomp(my $contents = `ls *\\@${runName}\\@fCounts.txt`);
        my @files = split(/[\\n]+/, $contents);
        foreach my $file (@files){
        $i++;
        $file=~/(.*)\\@${runName}\\@fCounts\\.txt/;
        my $libname = $1; 
        $a[$i]=$libname;
        open IN, $file;
            $_=<IN>;
            while(<IN>){
                my @v=split; 
                $b{$v[0]}{$i}=$v[$tf{$expectedCount}];
                $c{$v[0]}=$v[5]; #length column
            }
            close IN;
        }
        my $outfile="$runName"."_featureCounts.tsv";
        open OUT, ">$outfile";
        if ($runName eq "transcript_id") {
            print OUT "transcript\tlength";
        } else {
            print OUT "gene\tlength";
        }
    
        for(my $j=1;$j<=$i;$j++) {
            print OUT "\t$a[$j]";
        }
        print OUT "\n";
    
        foreach my $key (keys %b) {
            print OUT "$key\t$c{$key}";
            for(my $j=1;$j<=$i;$j++){
                print OUT "\t$b{$key}{$j}";
            }
            print OUT "\n";
        }
        close OUT;
         
    }
}


	

# Step 2: Merge summary files
for($l = 0; $l <= $#run_name; $l++) {
    my $runName = $run_name[$l];
    my @a=();
    my %b=();
    my $i=0;
    chomp(my $contents = `ls *\\@${runName}\\@fCounts.txt.summary`);
    my @files = split(/[\\n]+/, $contents);
    foreach my $file (@files){
        $i++;
        $file=~/(.*)\\@${runName}\\@fCounts\\.txt\\.summary/;
        my $libname = $1; 
        $a[$i]=$libname;
        open IN, $file;
        $_=<IN>;
        while(<IN>){
            my @v=split; 
            $b{$v[0]}{$i}=$v[1];
        }
        close IN;
    }
    my $outfile="$runName"."_featureCounts.sum.tsv";
    open OUT, ">$outfile";
    print OUT "criteria";
    for(my $j=1;$j<=$i;$j++) {
        print OUT "\t$a[$j]";
    }
    print OUT "\n";
    
    foreach my $key (keys %b) {
        print OUT "$key";
        for(my $j=1;$j<=$i;$j++){
            print OUT "\t$b{$key}{$j}";
        }
        print OUT "\n";
    }
    close OUT;
}

'''
}


//* autofill
//* platform
//* platform
//* autofill

process STAR_Module_merge_transcriptome_bam {

input:
 tuple val(oldname), file(bamfiles)

output:
 tuple val(oldname), file("${oldname}.bam")  ,emit:g14_15_merged_bams00 
 tuple val(oldname), file("*_sorted*bai")  ,emit:g14_15_bam_index11 
 tuple val(oldname), file("*_sorted*bam")  ,emit:g14_15_sorted_bam22 

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
    samtools merge !{oldname}.bam !{bamfiles.join(" ")} && samtools sort -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
else
    mv !{bamfiles.join(" ")} !{oldname}.bam 2>/dev/null || true
    samtools sort  -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
fi
'''
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 5
    $MEMORY = 30 
}
//* platform
//* platform
//* autofill

process DE_module_Prepare_DESeq2 {

input:
 path counts
 path groups_file
 path compare_file
 val run_DESeq2

output:
 path "DE_reports"  ,emit:g37_24_outputFile00_g37_37 
 val "_des"  ,emit:g37_24_postfix10_g37_33 
 path "DE_reports/outputs/*_all_deseq2_results.tsv"  ,emit:g37_24_outputFile21_g37_33 

container 'quay.io/viascientific/de_module:4.0'

when:
run_DESeq2 == 'yes'

script:

feature_type = params.DE_module_Prepare_DESeq2.feature_type

countsFile = ""
listOfFiles = counts.toString().split(" ")
if (listOfFiles.size() > 1){
	countsFile = listOfFiles[0]
	for(item in listOfFiles){
    	if (feature_type.equals('gene') && item.startsWith("gene") && item.toLowerCase().indexOf("tpm") < 0) {
    		countsFile = item
    		break;
    	} else if (feature_type.equals('transcript') && (item.startsWith("isoforms") || item.startsWith("transcript")) && item.toLowerCase().indexOf("tpm") < 0) {
    		countsFile = item
    		break;	
    	}
	}
} else {
	countsFile = counts
}

if (params.genome_build.startsWith("mousetest_")) {
	orgDB = 'org.Mm.eg.db'
} else if (params.genome_build.startsWith("human_")) {
	orgDB = 'org.Hs.eg.db'
} else if (params.genome_build.startsWith("mouse_")) {
	orgDB = 'org.Mm.eg.db'
} else if (params.genome_build.startsWith("rat_")) {
	orgDB = 'org.Rn.eg.db'
} else if (params.genome_build.startsWith("zebrafish_")) {
	orgDB = 'org.Dr.eg.db'
} else if (params.genome_build.startsWith("c_elegans_")) {
	orgDB = 'org.Ce.eg.db'
} else if (params.genome_build.startsWith("d_melanogaster_")) {
	orgDB = 'org.Dm.eg.db'
} else if (params.genome_build.startsWith("s_cerevisiae_")) {
	orgDB = 'org.Sc.sgd.db'
} else if (params.genome_build.startsWith("s_pombe_")) {
	orgDB = 'NA'
} else if (params.genome_build.startsWith("e_coli_")) {
	orgDB = 'org.EcK12.eg.db'
} else if (params.genome_build.startsWith("dog_")) {
	orgDB = 'org.Cf.eg.db'
} else if (params.genome_build == "custom"){
	orgDB = 'NA'
} else {
	orgDB = 'NA'
}

include_distribution = params.DE_module_Prepare_DESeq2.include_distribution
include_all2all = params.DE_module_Prepare_DESeq2.include_all2all
include_pca = params.DE_module_Prepare_DESeq2.include_pca

filter_type = params.DE_module_Prepare_DESeq2.filter_type
min_count = params.DE_module_Prepare_DESeq2.min_count
min_samples = params.DE_module_Prepare_DESeq2.min_samples
min_counts_per_sample = params.DE_module_Prepare_DESeq2.min_counts_per_sample
excluded_events = params.DE_module_Prepare_DESeq2.excluded_events

include_batch_correction = params.DE_module_Prepare_DESeq2.include_batch_correction
batch_correction_column = params.DE_module_Prepare_DESeq2.batch_correction_column
batch_correction_group_column = params.DE_module_Prepare_DESeq2.batch_correction_group_column
batch_normalization_algorithm = params.DE_module_Prepare_DESeq2.batch_normalization_algorithm

transformation = params.DE_module_Prepare_DESeq2.transformation
pca_color = params.DE_module_Prepare_DESeq2.pca_color
pca_shape = params.DE_module_Prepare_DESeq2.pca_shape
pca_fill = params.DE_module_Prepare_DESeq2.pca_fill
pca_transparency = params.DE_module_Prepare_DESeq2.pca_transparency
pca_label = params.DE_module_Prepare_DESeq2.pca_label

include_deseq2 = params.DE_module_Prepare_DESeq2.include_deseq2
input_mode = params.DE_module_Prepare_DESeq2.input_mode
design = params.DE_module_Prepare_DESeq2.design
fitType = params.DE_module_Prepare_DESeq2.fitType
use_batch_corrected_in_DE = params.DE_module_Prepare_DESeq2.use_batch_corrected_in_DE
apply_shrinkage = params.DE_module_Prepare_DESeq2.apply_shrinkage
shrinkage_type = params.DE_module_Prepare_DESeq2.shrinkage_type
include_volcano = params.DE_module_Prepare_DESeq2.include_volcano
include_ma = params.DE_module_Prepare_DESeq2.include_ma
include_heatmap = params.DE_module_Prepare_DESeq2.include_heatmap

padj_significance_cutoff = params.DE_module_Prepare_DESeq2.padj_significance_cutoff
fc_significance_cutoff = params.DE_module_Prepare_DESeq2.fc_significance_cutoff
padj_floor = params.DE_module_Prepare_DESeq2.padj_floor
fc_ceiling = params.DE_module_Prepare_DESeq2.fc_ceiling

convert_names = params.DE_module_Prepare_DESeq2.convert_names
count_file_names = params.DE_module_Prepare_DESeq2.count_file_names
converted_name = params.DE_module_Prepare_DESeq2.converted_name
num_labeled = params.DE_module_Prepare_DESeq2.num_labeled
highlighted_genes = params.DE_module_Prepare_DESeq2.highlighted_genes
include_volcano_highlighted = params.DE_module_Prepare_DESeq2.include_volcano_highlighted
include_ma_highlighted = params.DE_module_Prepare_DESeq2.include_ma_highlighted

//* @style @condition:{convert_names="true", count_file_names, converted_name},{convert_names="false"},{include_batch_correction="true", batch_correction_column, batch_correction_group_column, batch_normalization_algorithm, use_batch_corrected_in_DE},{include_batch_correction="false"},{include_deseq2="true", design, fitType, apply_shrinkage, shrinkage_type, include_volcano, include_ma, include_heatmap, padj_significance_cutoff, fc_significance_cutoff, padj_floor, fc_ceiling, convert_names, count_file_names, converted_name, num_labeled, highlighted_genes},{include_deseq2="false"},{apply_shrinkage="true", shrinkage_type},{apply_shrinkage="false"},{include_pca="true", transformation, pca_color, pca_shape, pca_fill, pca_transparency, pca_label},{include_pca="false"} @multicolumn:{include_distribution, include_all2all, include_pca},{filter_type, min_count, min_samples, min_counts_per_sample},{include_batch_correction, batch_correction_column, batch_correction_group_column, batch_normalization_algorithm},{design, fitType, use_batch_corrected_in_DE, apply_shrinkage, shrinkage_type},{pca_color, pca_shape, pca_fill, pca_transparency, pca_label},{padj_significance_cutoff, fc_significance_cutoff, padj_floor, fc_ceiling},{convert_names, count_file_names, converted_name},{include_volcano, include_ma, include_heatmap}{include_volcano_highlighted,include_ma_highlighted} 

include_distribution = include_distribution == 'true' ? 'TRUE' : 'FALSE'
include_all2all = include_all2all == 'true' ? 'TRUE' : 'FALSE'
include_pca = include_pca == 'true' ? 'TRUE' : 'FALSE'
include_batch_correction = include_batch_correction == 'true' ? 'TRUE' : 'FALSE'
include_deseq2 = include_deseq2 == 'true' ? 'TRUE' : 'FALSE'
use_batch_corrected_in_DE = use_batch_corrected_in_DE  == 'true' ? 'TRUE' : 'FALSE'
apply_shrinkage = apply_shrinkage == 'true' ? 'TRUE' : 'FALSE'
convert_names = convert_names == 'true' ? 'TRUE' : 'FALSE'
include_ma = include_ma == 'true' ? 'TRUE' : 'FALSE'
include_volcano = include_volcano == 'true' ? 'TRUE' : 'FALSE'
include_heatmap = include_heatmap == 'true' ? 'TRUE' : 'FALSE'

excluded_events = excluded_events.replace("\n", " ").replace(',', ' ')

pca_color_arg = pca_color.equals('') ? '' : '--pca-color ' + pca_color
pca_shape_arg = pca_shape.equals('') ? '' : '--pca-shape ' + pca_shape
pca_fill_arg = pca_fill.equals('') ? '' : '--pca-fill ' + pca_fill
pca_transparency_arg = pca_transparency.equals('') ? '' : '--pca-transparency ' + pca_transparency
pca_label_arg = pca_label.equals('') ? '' : '--pca-label ' + pca_label

highlighted_genes = highlighted_genes.replace("\n", " ").replace(',', ' ')

excluded_events_arg = excluded_events.equals('') ? '' : '--excluded-events ' + excluded_events
highlighted_genes_arg = highlighted_genes.equals('') ? '' : '--highlighted-genes ' + highlighted_genes
include_ma_highlighted = include_ma_highlighted == 'true' ? 'TRUE' : 'FALSE'
include_volcano_highlighted = include_volcano_highlighted == 'true' ? 'TRUE' : 'FALSE'

threads = task.cpus

"""
mkdir reports
mkdir inputs
mkdir outputs
cp ${groups_file} inputs/${groups_file}
cp ${compare_file} inputs/${compare_file}
cp ${countsFile} inputs/${countsFile}
cp inputs/${groups_file} inputs/.de_metadata.txt
cp inputs/${compare_file} inputs/.comparisons.tsv

prepare_DESeq2.py \
--counts inputs/${countsFile} --groups inputs/${groups_file} --comparisons inputs/${compare_file} --feature-type ${feature_type} \
--include-distribution ${include_distribution} --include-all2all ${include_all2all} --include-pca ${include_pca} \
--filter-type ${filter_type} --min-counts-per-event ${min_count} --min-samples-per-event ${min_samples} --min-counts-per-sample ${min_counts_per_sample} \
--transformation ${transformation} ${pca_color_arg} ${pca_shape_arg} ${pca_fill_arg} ${pca_transparency_arg} ${pca_label_arg} \
--include-batch-correction ${include_batch_correction} --batch-correction-column ${batch_correction_column} --batch-correction-group-column ${batch_correction_group_column} --batch-normalization-algorithm ${batch_normalization_algorithm} \
--include-DESeq2 ${include_deseq2} --input-mode ${input_mode} --design '${design}' --fitType ${fitType} --use-batch-correction-in-DE ${use_batch_corrected_in_DE} --apply-shrinkage ${apply_shrinkage} --shrinkage-type ${shrinkage_type} \
--include-volcano ${include_volcano} --include-ma ${include_ma} --include-heatmap ${include_heatmap} \
--padj-significance-cutoff ${padj_significance_cutoff} --fc-significance-cutoff ${fc_significance_cutoff} --padj-floor ${padj_floor} --fc-ceiling ${fc_ceiling} \
--convert-names ${convert_names} --count-file-names ${count_file_names} --converted-names ${converted_name} --org-db ${orgDB} --num-labeled ${num_labeled} \
${highlighted_genes_arg} --include-volcano-highlighted ${include_volcano_highlighted} --include-ma-highlighted ${include_ma_highlighted} \
${excluded_events_arg} \
--threads ${threads}

mkdir DE_reports
mv *.Rmd DE_reports
mv *.html DE_reports
mv inputs DE_reports/
mv outputs DE_reports/
"""

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 5
    $MEMORY = 30 
}
//* platform
//* platform
//* autofill

process DE_module_Prepare_GSEA_DESeq2 {

input:
 val postfix
 path input
 val run_GSEA

output:
 path "GSEA_reports"  ,emit:g37_33_outputFile01_g37_37 
 path "GSEA"  ,emit:g37_33_outputFile11 

container 'quay.io/viascientific/gsea_module:2.0.0'

// SET second output to "{GSEA,outputs}" when launched apps can reach parent directory

when:
run_GSEA == 'yes'

script:

if (postfix.equals(" ")) {
	postfix = ''
}

event_column = params.DE_module_Prepare_GSEA_DESeq2.event_column
fold_change_column = params.DE_module_Prepare_GSEA_DESeq2.fold_change_column

if (params.genome_build.startsWith("human_")) {
	species = 'human'
} else if (params.genome_build.startsWith("mouse_")) {
	species = 'mouse'
} else {
	species = 'NA'
}

local_species = params.DE_module_Prepare_GSEA_DESeq2.local_species
H = params.DE_module_Prepare_GSEA_DESeq2.H
C1 = params.DE_module_Prepare_GSEA_DESeq2.C1
C2 = params.DE_module_Prepare_GSEA_DESeq2.C2
C2_CGP = params.DE_module_Prepare_GSEA_DESeq2.C2_CGP
C2_CP = params.DE_module_Prepare_GSEA_DESeq2.C2_CP
C3_MIR = params.DE_module_Prepare_GSEA_DESeq2.C3_MIR
C3_TFT = params.DE_module_Prepare_GSEA_DESeq2.C3_TFT
C4 = params.DE_module_Prepare_GSEA_DESeq2.C4
C4_CGN = params.DE_module_Prepare_GSEA_DESeq2.C4_CGN
C4_CM = params.DE_module_Prepare_GSEA_DESeq2.C4_CM
C5 = params.DE_module_Prepare_GSEA_DESeq2.C5
C5_GO = params.DE_module_Prepare_GSEA_DESeq2.C5_GO
C5_GO_BP = params.DE_module_Prepare_GSEA_DESeq2.C5_GO_BP
C5_GO_CC = params.DE_module_Prepare_GSEA_DESeq2.C5_GO_CC
C5_GO_MF = params.DE_module_Prepare_GSEA_DESeq2.C5_GO_MF
C5_HPO = params.DE_module_Prepare_GSEA_DESeq2.C5_HPO
C6 = params.DE_module_Prepare_GSEA_DESeq2.C6
C7 = params.DE_module_Prepare_GSEA_DESeq2.C7
C8 = params.DE_module_Prepare_GSEA_DESeq2.C8
MH = params.DE_module_Prepare_GSEA_DESeq2.MH
M1 = params.DE_module_Prepare_GSEA_DESeq2.M1
M2 = params.DE_module_Prepare_GSEA_DESeq2.M2
M2_CGP = params.DE_module_Prepare_GSEA_DESeq2.M2_CGP
M2_CP = params.DE_module_Prepare_GSEA_DESeq2.M2_CP
M3_GTRD = params.DE_module_Prepare_GSEA_DESeq2.M3_GTRD
M3_miRDB = params.DE_module_Prepare_GSEA_DESeq2.M3_miRDB
M5 = params.DE_module_Prepare_GSEA_DESeq2.M5
M5_GO = params.DE_module_Prepare_GSEA_DESeq2.M5_GO
M5_GO_BP = params.DE_module_Prepare_GSEA_DESeq2.M5_GO_BP
M5_GO_CC = params.DE_module_Prepare_GSEA_DESeq2.M5_GO_CC
M5_GO_MF = params.DE_module_Prepare_GSEA_DESeq2.M5_GO_MF
M5_MPT = params.DE_module_Prepare_GSEA_DESeq2.M5_MPT
M8 = params.DE_module_Prepare_GSEA_DESeq2.M8

minSize = params.DE_module_Prepare_GSEA_DESeq2.minSize
maxSize = params.DE_module_Prepare_GSEA_DESeq2.maxSize

nes_significance_cutoff = params.DE_module_Prepare_GSEA_DESeq2.nes_significance_cutoff
padj_significance_cutoff = params.DE_module_Prepare_GSEA_DESeq2.padj_significance_cutoff

seed = params.DE_module_Prepare_GSEA_DESeq2.seed

//* @style @condition:{local_species="human",H,C1,C2,C2_CGP,C2_CP,C3_MIR,C3_TFT,C4,C4_CGN,C4_CM,C5,C5_GO,C5_GO_BP,C5_GO_CC,C5_GO_MF,C5_HPO,C6,C7,C8},{local_species="mouse",MH,M1,M2,M2_CGP,M2_CP,M3_GTRD,M3_miRDB,M5,M5_GO,M5_GO_BP,M5_GO_CC,M5_GO_MF,M5_MPT,M8} @multicolumn:{event_column, fold_change_column}, {gmt_list, minSize, maxSize},{H,C1,C2,C2_CGP}, {C2_CP,C3_MIR,C3_TFT,C4},{C4_CGN,C4_CM, C5,C5_GO},{C5_GO_BP,C5_GO_CC,C5_GO_MF,C5_HPO},{C6,C7,C8},{MH,M1,M2,M2_CGP},{M2_CP,M3_GTRD,M3_miRDB,M5},{M5_GO,M5_GO_BP,M5_GO_CC,M5_GO_MF},{M5_MPT,M8},{nes_significance_cutoff, padj_significance_cutoff}

H        = H        == 'true' && species == 'human' ? ' h.all.v2023.2.Hs.symbols.gmt'    : ''
C1       = C1       == 'true' && species == 'human' ? ' c1.all.v2023.2.Hs.symbols.gmt'   : ''
C2       = C2       == 'true' && species == 'human' ? ' c2.all.v2023.2.Hs.symbols.gmt'   : ''
C2_CGP   = C2_CGP   == 'true' && species == 'human' ? ' c2.cgp.v2023.2.Hs.symbols.gmt'   : ''
C2_CP    = C2_CP    == 'true' && species == 'human' ? ' c2.cp.v2023.2.Hs.symbols.gmt'    : ''
C3_MIR   = C3_MIR   == 'true' && species == 'human' ? ' c3.mir.v2023.2.Hs.symbols.gmt'   : ''
C3_TFT   = C3_TFT   == 'true' && species == 'human' ? ' c3.tft.v2023.2.Hs.symbols.gmt'   : ''
C4       = C4       == 'true' && species == 'human' ? ' c4.all.v2023.2.Hs.symbols.gmt'   : ''
C4_CGN   = C4_CGN   == 'true' && species == 'human' ? ' c4.cgn.v2023.2.Hs.symbols.gmt'   : ''
C4_CM    = C4_CM    == 'true' && species == 'human' ? ' c4.cm.v2023.2.Hs.symbols.gmt'    : ''
C5       = C5       == 'true' && species == 'human' ? ' c5.all.v2023.2.Hs.symbols.gmt'   : ''
C5_GO    = C5_GO    == 'true' && species == 'human' ? ' c5.go.v2023.2.Hs.symbols.gmt'    : ''
C5_GO_BP = C5_GO_BP == 'true' && species == 'human' ? ' c5.go.bp.v2023.2.Hs.symbols.gmt' : ''
C5_GO_CC = C5_GO_CC == 'true' && species == 'human' ? ' c5.go.cc.v2023.2.Hs.symbols.gmt' : ''
C5_GO_MF = C5_GO_MF == 'true' && species == 'human' ? ' c5.go.mf.v2023.2.Hs.symbols.gmt' : ''
C5_HPO   = C5_HPO   == 'true' && species == 'human' ? ' c5.hpo.v2023.2.Hs.symbols.gmt'   : ''
C6       = C6       == 'true' && species == 'human' ? ' c6.all.v2023.2.Hs.symbols.gmt'   : ''
C7       = C7       == 'true' && species == 'human' ? ' c7.all.v2023.2.Hs.symbols.gmt'   : ''
C8       = C8       == 'true' && species == 'human' ? ' c8.all.v2023.2.Hs.symbols.gmt'   : ''
MH       = MH       == 'true' && species == 'mouse' ? ' mh.all.v2023.2.Mm.symbols.gmt'   : ''    
M1       = M1       == 'true' && species == 'mouse' ? ' m1.all.v2023.2.Mm.symbols.gmt'   : ''    
M2       = M2       == 'true' && species == 'mouse' ? ' m2.all.v2023.2.Mm.symbols.gmt'   : ''    
M2_CGP   = M2_CGP   == 'true' && species == 'mouse' ? ' m2.cgp.v2023.2.Mm.symbols.gmt'   : ''
M2_CP    = M2_CP    == 'true' && species == 'mouse' ? ' m2.cp.v2023.2.Mm.symbols.gmt'    : '' 
M3_GTRD  = M3_GTRD  == 'true' && species == 'mouse' ? ' m3.gtrd.v2023.2.Mm.symbols.gmt'  : ''
M3_miRDB = M3_miRDB == 'true' && species == 'mouse' ? ' m3.mirdb.v2023.2.Mm.symbols.gmt' : ''
M5       = M5       == 'true' && species == 'mouse' ? ' m5.all.v2023.2.Mm.symbols.gmt'   : ''    
M5_GO    = M5_GO    == 'true' && species == 'mouse' ? ' m5.go.v2023.2.Mm.symbols.gmt'    : '' 
M5_GO_BP = M5_GO_BP == 'true' && species == 'mouse' ? ' m5.go.bp.v2023.2.Mm.symbols.gmt' : ''
M5_GO_CC = M5_GO_CC == 'true' && species == 'mouse' ? ' m5.go.cc.v2023.2.Mm.symbols.gmt' : ''
M5_GO_MF = M5_GO_MF == 'true' && species == 'mouse' ? ' m5.go.mf.v2023.2.Mm.symbols.gmt' : ''
M5_MPT   = M5_MPT   == 'true' && species == 'mouse' ? ' m5.mpt.v2023.2.Mm.symbols.gmt'   : ''
M8       = M8       == 'true' && species == 'mouse' ? ' m8.all.v2023.2.Mm.symbols.gmt'   : ''

gmt_list = H + C1 + C2 + C2_CGP + C2_CP + C3_MIR + C3_TFT + C4 + C4_CGN + C4_CM + C5 + C5_GO + C5_GO_BP + C5_GO_CC + C5_GO_MF + C5_HPO + C6 + C7 + C8 + MH + M1 + M2 + M2_CGP + M2_CP + M3_GTRD + M3_miRDB + M5 + M5_GO + M5_GO_BP + M5_GO_CC + M5_GO_MF + M5_MPT + M8

if (gmt_list.equals("")){
	gmt_list = 'h.all.v2023.2.Hs.symbols.gmt'
}

"""
prepare_GSEA.py \
--input ${input} --species ${species} --event-column ${event_column} --fold-change-column ${fold_change_column} \
--GMT-key /data/gmt_key.txt --GMT-source /data --GMT-list ${gmt_list} --minSize ${minSize} --maxSize ${maxSize} \
--NES ${nes_significance_cutoff} --pvalue ${padj_significance_cutoff} \
--seed ${seed} --threads ${task.cpus} --postfix '_gsea${postfix}'

cp -R outputs GSEA/

mkdir GSEA_reports
cp -R GSEA GSEA_reports/
"""

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 4 
}
//* platform
//* platform
//* autofill

process DE_module_DESeq2_Analysis {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /(.*.Rmd|.*.html|inputs|outputs|GSEA)$/) "DESeq2/$filename"}
input:
 path DE_reports
 path GSEA_reports

output:
 path "{*.Rmd,*.html,inputs,outputs,GSEA}"  ,emit:g37_37_outputDir00 

container 'quay.io/viascientific/de_module:2.1.0'

script:

"""
mv DE_reports/* .

if [ -d "GSEA_reports" ]; then
    mv GSEA_reports/* .
fi
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 5
    $MEMORY = 30 
}
//* platform
//* platform
//* autofill

process DE_module_Prepare_LimmaVoom {

input:
 path counts
 path groups_file
 path compare_file
 val run_limmaVoom

output:
 path "DE_reports"  ,emit:g37_42_outputFile00_g37_39 
 val "_lv"  ,emit:g37_42_postfix10_g37_41 
 path "DE_reports/outputs/*_all_limmaVoom_results.tsv"  ,emit:g37_42_outputFile21_g37_41 

container 'quay.io/viascientific/de_module:4.0'

when:
run_limmaVoom == 'yes'

script:

feature_type = params.DE_module_Prepare_LimmaVoom.feature_type

countsFile = ""
listOfFiles = counts.toString().split(" ")
if (listOfFiles.size() > 1){
	countsFile = listOfFiles[0]
	for(item in listOfFiles){
    	if (feature_type.equals('gene') && item.startsWith("gene") && item.toLowerCase().indexOf("tpm") < 0) {
    		countsFile = item
    		break;
    	} else if (feature_type.equals('transcript') && (item.startsWith("isoforms") || item.startsWith("transcript")) && item.toLowerCase().indexOf("tpm") < 0) {
    		countsFile = item
    		break;	
    	}
	}
} else {
	countsFile = counts
}

if (params.genome_build.startsWith("mousetest_")) {
	orgDB = 'org.Mm.eg.db'
} else if (params.genome_build.startsWith("human_")) {
	orgDB = 'org.Hs.eg.db'
} else if (params.genome_build.startsWith("mouse_")) {
	orgDB = 'org.Mm.eg.db'
} else if (params.genome_build.startsWith("rat_")) {
	orgDB = 'org.Rn.eg.db'
} else if (params.genome_build.startsWith("zebrafish_")) {
	orgDB = 'org.Dr.eg.db'
} else if (params.genome_build.startsWith("c_elegans_")) {
	orgDB = 'org.Ce.eg.db'
} else if (params.genome_build.startsWith("d_melanogaster_")) {
	orgDB = 'org.Dm.eg.db'
} else if (params.genome_build.startsWith("s_cerevisiae_")) {
	orgDB = 'org.Sc.sgd.db'
} else if (params.genome_build.startsWith("s_pombe_")) {
	orgDB = 'NA'
} else if (params.genome_build.startsWith("e_coli_")) {
	orgDB = 'org.EcK12.eg.db'
} else if (params.genome_build.startsWith("dog_")) {
	orgDB = 'org.Cf.eg.db'
} else if (params.genome_build == "custom"){
	orgDB = 'NA'
} else {
	orgDB = 'NA'
}

include_distribution = params.DE_module_Prepare_LimmaVoom.include_distribution
include_all2all = params.DE_module_Prepare_LimmaVoom.include_all2all
include_pca = params.DE_module_Prepare_LimmaVoom.include_pca

filter_type = params.DE_module_Prepare_LimmaVoom.filter_type
min_count = params.DE_module_Prepare_LimmaVoom.min_count
min_samples = params.DE_module_Prepare_LimmaVoom.min_samples
min_counts_per_sample = params.DE_module_Prepare_LimmaVoom.min_counts_per_sample
excluded_events = params.DE_module_Prepare_LimmaVoom.excluded_events

include_batch_correction = params.DE_module_Prepare_LimmaVoom.include_batch_correction
batch_correction_column = params.DE_module_Prepare_LimmaVoom.batch_correction_column
batch_correction_group_column = params.DE_module_Prepare_LimmaVoom.batch_correction_group_column
batch_normalization_algorithm = params.DE_module_Prepare_LimmaVoom.batch_normalization_algorithm

transformation = params.DE_module_Prepare_LimmaVoom.transformation
pca_color = params.DE_module_Prepare_LimmaVoom.pca_color
pca_shape = params.DE_module_Prepare_LimmaVoom.pca_shape
pca_fill = params.DE_module_Prepare_LimmaVoom.pca_fill
pca_transparency = params.DE_module_Prepare_LimmaVoom.pca_transparency
pca_label = params.DE_module_Prepare_LimmaVoom.pca_label

include_limma = params.DE_module_Prepare_LimmaVoom.include_limma
use_batch_corrected_in_DE = params.DE_module_Prepare_LimmaVoom.use_batch_corrected_in_DE
normalization_method = params.DE_module_Prepare_LimmaVoom.normalization_method
logratioTrim = params.DE_module_Prepare_LimmaVoom.logratioTrim
sumTrim = params.DE_module_Prepare_LimmaVoom.sumTrim
Acutoff = params.DE_module_Prepare_LimmaVoom.Acutoff
doWeighting = params.DE_module_Prepare_LimmaVoom.doWeighting
p = params.DE_module_Prepare_LimmaVoom.p
include_volcano = params.DE_module_Prepare_LimmaVoom.include_volcano
include_ma = params.DE_module_Prepare_LimmaVoom.include_ma
include_heatmap = params.DE_module_Prepare_LimmaVoom.include_heatmap

padj_significance_cutoff = params.DE_module_Prepare_LimmaVoom.padj_significance_cutoff
fc_significance_cutoff = params.DE_module_Prepare_LimmaVoom.fc_significance_cutoff
padj_floor = params.DE_module_Prepare_LimmaVoom.padj_floor
fc_ceiling = params.DE_module_Prepare_LimmaVoom.fc_ceiling

convert_names = params.DE_module_Prepare_LimmaVoom.convert_names
count_file_names = params.DE_module_Prepare_LimmaVoom.count_file_names
converted_name = params.DE_module_Prepare_LimmaVoom.converted_name
num_labeled = params.DE_module_Prepare_LimmaVoom.num_labeled
highlighted_genes = params.DE_module_Prepare_LimmaVoom.highlighted_genes
include_volcano_highlighted = params.DE_module_Prepare_LimmaVoom.include_volcano_highlighted
include_ma_highlighted = params.DE_module_Prepare_LimmaVoom.include_ma_highlighted

//* @style @condition:{convert_names="true", count_file_names, converted_name},{convert_names="false"},{include_batch_correction="true", batch_correction_column, batch_correction_group_column, batch_normalization_algorithm,use_batch_corrected_in_DE},{include_batch_correction="false"},{include_limma="true", normalization_method, logratioTrim, sumTrim, doWeighting, Acutoff, include_volcano, include_ma, include_heatmap, padj_significance_cutoff, fc_significance_cutoff, padj_floor, fc_ceiling, convert_names, count_file_names, converted_name, num_labeled, highlighted_genes},{include_limma="false"},{normalization_method="TMM", logratioTrim, sumTrim, doWeighting, Acutoff},{normalization_method="TMMwsp", logratioTrim, sumTrim, doWeighting, Acutoff},{normalization_method="RLE"},{normalization_method="upperquartile", p},{normalization_method="none"},{include_pca="true", transformation, pca_color, pca_shape, pca_fill, pca_transparency, pca_label},{include_pca="false"} @multicolumn:{include_distribution, include_all2all, include_pca},{filter_type, min_count, min_samples,min_counts_per_sample},{include_batch_correction, batch_correction_column, batch_correction_group_column, batch_normalization_algorithm},{pca_color, pca_shape, pca_fill, pca_transparency, pca_label},{include_limma, use_batch_corrected_in_DE},{normalization_method,logratioTrim,sumTrim,doWeighting,Acutoff,p},{padj_significance_cutoff, fc_significance_cutoff, padj_floor, fc_ceiling},{convert_names, count_file_names, converted_name},{include_volcano, include_ma, include_heatmap}{include_volcano_highlighted,include_ma_highlighted} 

include_distribution = include_distribution == 'true' ? 'TRUE' : 'FALSE'
include_all2all = include_all2all == 'true' ? 'TRUE' : 'FALSE'
include_pca = include_pca == 'true' ? 'TRUE' : 'FALSE'
include_batch_correction = include_batch_correction == 'true' ? 'TRUE' : 'FALSE'
include_limma = include_limma == 'true' ? 'TRUE' : 'FALSE'
use_batch_corrected_in_DE = use_batch_corrected_in_DE  == 'true' ? 'TRUE' : 'FALSE'
convert_names = convert_names == 'true' ? 'TRUE' : 'FALSE'
include_ma = include_ma == 'true' ? 'TRUE' : 'FALSE'
include_volcano = include_volcano == 'true' ? 'TRUE' : 'FALSE'
include_heatmap = include_heatmap == 'true' ? 'TRUE' : 'FALSE'

doWeighting = doWeighting == 'true' ? 'TRUE' : 'FALSE'
TMM_args = normalization_method.equals('TMM') || normalization_method.equals('TMMwsp') ? '--logratio-trim ' + logratioTrim + ' --sum-trim ' + sumTrim + ' --do-weighting ' + doWeighting + ' --A-cutoff="' + Acutoff + '"' : ''
upperquartile_args = normalization_method.equals('upperquartile') ? '--p ' + p : ''

excluded_events = excluded_events.replace("\n", " ").replace(',', ' ')

pca_color_arg = pca_color.equals('') ? '' : '--pca-color ' + pca_color
pca_shape_arg = pca_shape.equals('') ? '' : '--pca-shape ' + pca_shape
pca_fill_arg = pca_fill.equals('') ? '' : '--pca-fill ' + pca_fill
pca_transparency_arg = pca_transparency.equals('') ? '' : '--pca-transparency ' + pca_transparency
pca_label_arg = pca_label.equals('') ? '' : '--pca-label ' + pca_label

highlighted_genes = highlighted_genes.replace("\n", " ").replace(',', ' ')

excluded_events_arg = excluded_events.equals('') ? '' : '--excluded-events ' + excluded_events
highlighted_genes_arg = highlighted_genes.equals('') ? '' : '--highlighted-genes ' + highlighted_genes
include_ma_highlighted = include_ma_highlighted == 'true' ? 'TRUE' : 'FALSE'
include_volcano_highlighted = include_volcano_highlighted == 'true' ? 'TRUE' : 'FALSE'

threads = task.cpus

"""
mkdir inputs
mkdir outputs
cp ${groups_file} inputs/${groups_file}
cp ${compare_file} inputs/${compare_file}
cp ${countsFile} inputs/${countsFile}
cp inputs/${groups_file} inputs/.de_metadata.txt
cp inputs/${compare_file} inputs/.comparisons.tsv

prepare_limmaVoom.py \
--counts inputs/${countsFile} --groups inputs/${groups_file} --comparisons inputs/${compare_file} --feature-type ${feature_type} \
--include-distribution ${include_distribution} --include-all2all ${include_all2all} --include-pca ${include_pca} \
--filter-type ${filter_type} --min-counts-per-event ${min_count} --min-samples-per-event ${min_samples} --min-counts-per-sample ${min_counts_per_sample} \
--transformation ${transformation} ${pca_color_arg} ${pca_shape_arg} ${pca_fill_arg} ${pca_transparency_arg} ${pca_label_arg} \
--include-batch-correction ${include_batch_correction} --batch-correction-column ${batch_correction_column} --batch-correction-group-column ${batch_correction_group_column} --batch-normalization-algorithm ${batch_normalization_algorithm} \
--include-limma ${include_limma} \
--use-batch-correction-in-DE ${use_batch_corrected_in_DE} --normalization-method ${normalization_method} ${TMM_args} ${upperquartile_args} \
--include-volcano ${include_volcano} --include-ma ${include_ma} --include-heatmap ${include_heatmap} \
--padj-significance-cutoff ${padj_significance_cutoff} --fc-significance-cutoff ${fc_significance_cutoff} --padj-floor ${padj_floor} --fc-ceiling ${fc_ceiling} \
--convert-names ${convert_names} --count-file-names ${count_file_names} --org-db ${orgDB} --converted-names ${converted_name} --num-labeled ${num_labeled} \
${highlighted_genes_arg} --include-volcano-highlighted ${include_volcano_highlighted} --include-ma-highlighted ${include_ma_highlighted} \
${excluded_events_arg} \
--threads ${threads}

mkdir DE_reports
mv *.Rmd DE_reports
mv *.html DE_reports
mv inputs DE_reports/
mv outputs DE_reports/
"""

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 5
    $MEMORY = 30 
}
//* platform
//* platform
//* autofill

process DE_module_Prepare_GSEA_LimmaVoom {

input:
 val postfix
 path input
 val run_GSEA

output:
 path "GSEA_reports"  ,emit:g37_41_outputFile01_g37_39 
 path "GSEA"  ,emit:g37_41_outputFile11 

container 'quay.io/viascientific/gsea_module:2.0.0'

// SET second output to "{GSEA,outputs}" when launched apps can reach parent directory

when:
run_GSEA == 'yes'

script:

if (postfix.equals(" ")) {
	postfix = ''
}

event_column = params.DE_module_Prepare_GSEA_LimmaVoom.event_column
fold_change_column = params.DE_module_Prepare_GSEA_LimmaVoom.fold_change_column

if (params.genome_build.startsWith("human_")) {
	species = 'human'
} else if (params.genome_build.startsWith("mouse_")) {
	species = 'mouse'
} else {
	species = 'NA'
}

local_species = params.DE_module_Prepare_GSEA_LimmaVoom.local_species
H = params.DE_module_Prepare_GSEA_LimmaVoom.H
C1 = params.DE_module_Prepare_GSEA_LimmaVoom.C1
C2 = params.DE_module_Prepare_GSEA_LimmaVoom.C2
C2_CGP = params.DE_module_Prepare_GSEA_LimmaVoom.C2_CGP
C2_CP = params.DE_module_Prepare_GSEA_LimmaVoom.C2_CP
C3_MIR = params.DE_module_Prepare_GSEA_LimmaVoom.C3_MIR
C3_TFT = params.DE_module_Prepare_GSEA_LimmaVoom.C3_TFT
C4 = params.DE_module_Prepare_GSEA_LimmaVoom.C4
C4_CGN = params.DE_module_Prepare_GSEA_LimmaVoom.C4_CGN
C4_CM = params.DE_module_Prepare_GSEA_LimmaVoom.C4_CM
C5 = params.DE_module_Prepare_GSEA_LimmaVoom.C5
C5_GO = params.DE_module_Prepare_GSEA_LimmaVoom.C5_GO
C5_GO_BP = params.DE_module_Prepare_GSEA_LimmaVoom.C5_GO_BP
C5_GO_CC = params.DE_module_Prepare_GSEA_LimmaVoom.C5_GO_CC
C5_GO_MF = params.DE_module_Prepare_GSEA_LimmaVoom.C5_GO_MF
C5_HPO = params.DE_module_Prepare_GSEA_LimmaVoom.C5_HPO
C6 = params.DE_module_Prepare_GSEA_LimmaVoom.C6
C7 = params.DE_module_Prepare_GSEA_LimmaVoom.C7
C8 = params.DE_module_Prepare_GSEA_LimmaVoom.C8
MH = params.DE_module_Prepare_GSEA_LimmaVoom.MH
M1 = params.DE_module_Prepare_GSEA_LimmaVoom.M1
M2 = params.DE_module_Prepare_GSEA_LimmaVoom.M2
M2_CGP = params.DE_module_Prepare_GSEA_LimmaVoom.M2_CGP
M2_CP = params.DE_module_Prepare_GSEA_LimmaVoom.M2_CP
M3_GTRD = params.DE_module_Prepare_GSEA_LimmaVoom.M3_GTRD
M3_miRDB = params.DE_module_Prepare_GSEA_LimmaVoom.M3_miRDB
M5 = params.DE_module_Prepare_GSEA_LimmaVoom.M5
M5_GO = params.DE_module_Prepare_GSEA_LimmaVoom.M5_GO
M5_GO_BP = params.DE_module_Prepare_GSEA_LimmaVoom.M5_GO_BP
M5_GO_CC = params.DE_module_Prepare_GSEA_LimmaVoom.M5_GO_CC
M5_GO_MF = params.DE_module_Prepare_GSEA_LimmaVoom.M5_GO_MF
M5_MPT = params.DE_module_Prepare_GSEA_LimmaVoom.M5_MPT
M8 = params.DE_module_Prepare_GSEA_LimmaVoom.M8

minSize = params.DE_module_Prepare_GSEA_LimmaVoom.minSize
maxSize = params.DE_module_Prepare_GSEA_LimmaVoom.maxSize

nes_significance_cutoff = params.DE_module_Prepare_GSEA_LimmaVoom.nes_significance_cutoff
padj_significance_cutoff = params.DE_module_Prepare_GSEA_LimmaVoom.padj_significance_cutoff

seed = params.DE_module_Prepare_GSEA_LimmaVoom.seed

//* @style @condition:{local_species="human",H,C1,C2,C2_CGP,C2_CP,C3_MIR,C3_TFT,C4,C4_CGN,C4_CM,C5,C5_GO,C5_GO_BP,C5_GO_CC,C5_GO_MF,C5_HPO,C6,C7,C8},{local_species="mouse",MH,M1,M2,M2_CGP,M2_CP,M3_GTRD,M3_miRDB,M5,M5_GO,M5_GO_BP,M5_GO_CC,M5_GO_MF,M5_MPT,M8} @multicolumn:{event_column, fold_change_column}, {gmt_list, minSize, maxSize},{H,C1,C2,C2_CGP}, {C2_CP,C3_MIR,C3_TFT,C4},{C4_CGN,C4_CM, C5,C5_GO},{C5_GO_BP,C5_GO_CC,C5_GO_MF,C5_HPO},{C6,C7,C8},{MH,M1,M2,M2_CGP},{M2_CP,M3_GTRD,M3_miRDB,M5},{M5_GO,M5_GO_BP,M5_GO_CC,M5_GO_MF},{M5_MPT,M8},{nes_significance_cutoff, padj_significance_cutoff}

H        = H        == 'true' && species == 'human' ? ' h.all.v2023.2.Hs.symbols.gmt'    : ''
C1       = C1       == 'true' && species == 'human' ? ' c1.all.v2023.2.Hs.symbols.gmt'   : ''
C2       = C2       == 'true' && species == 'human' ? ' c2.all.v2023.2.Hs.symbols.gmt'   : ''
C2_CGP   = C2_CGP   == 'true' && species == 'human' ? ' c2.cgp.v2023.2.Hs.symbols.gmt'   : ''
C2_CP    = C2_CP    == 'true' && species == 'human' ? ' c2.cp.v2023.2.Hs.symbols.gmt'    : ''
C3_MIR   = C3_MIR   == 'true' && species == 'human' ? ' c3.mir.v2023.2.Hs.symbols.gmt'   : ''
C3_TFT   = C3_TFT   == 'true' && species == 'human' ? ' c3.tft.v2023.2.Hs.symbols.gmt'   : ''
C4       = C4       == 'true' && species == 'human' ? ' c4.all.v2023.2.Hs.symbols.gmt'   : ''
C4_CGN   = C4_CGN   == 'true' && species == 'human' ? ' c4.cgn.v2023.2.Hs.symbols.gmt'   : ''
C4_CM    = C4_CM    == 'true' && species == 'human' ? ' c4.cm.v2023.2.Hs.symbols.gmt'    : ''
C5       = C5       == 'true' && species == 'human' ? ' c5.all.v2023.2.Hs.symbols.gmt'   : ''
C5_GO    = C5_GO    == 'true' && species == 'human' ? ' c5.go.v2023.2.Hs.symbols.gmt'    : ''
C5_GO_BP = C5_GO_BP == 'true' && species == 'human' ? ' c5.go.bp.v2023.2.Hs.symbols.gmt' : ''
C5_GO_CC = C5_GO_CC == 'true' && species == 'human' ? ' c5.go.cc.v2023.2.Hs.symbols.gmt' : ''
C5_GO_MF = C5_GO_MF == 'true' && species == 'human' ? ' c5.go.mf.v2023.2.Hs.symbols.gmt' : ''
C5_HPO   = C5_HPO   == 'true' && species == 'human' ? ' c5.hpo.v2023.2.Hs.symbols.gmt'   : ''
C6       = C6       == 'true' && species == 'human' ? ' c6.all.v2023.2.Hs.symbols.gmt'   : ''
C7       = C7       == 'true' && species == 'human' ? ' c7.all.v2023.2.Hs.symbols.gmt'   : ''
C8       = C8       == 'true' && species == 'human' ? ' c8.all.v2023.2.Hs.symbols.gmt'   : ''
MH       = MH       == 'true' && species == 'mouse' ? ' mh.all.v2023.2.Mm.symbols.gmt'   : ''    
M1       = M1       == 'true' && species == 'mouse' ? ' m1.all.v2023.2.Mm.symbols.gmt'   : ''    
M2       = M2       == 'true' && species == 'mouse' ? ' m2.all.v2023.2.Mm.symbols.gmt'   : ''    
M2_CGP   = M2_CGP   == 'true' && species == 'mouse' ? ' m2.cgp.v2023.2.Mm.symbols.gmt'   : ''
M2_CP    = M2_CP    == 'true' && species == 'mouse' ? ' m2.cp.v2023.2.Mm.symbols.gmt'    : '' 
M3_GTRD  = M3_GTRD  == 'true' && species == 'mouse' ? ' m3.gtrd.v2023.2.Mm.symbols.gmt'  : ''
M3_miRDB = M3_miRDB == 'true' && species == 'mouse' ? ' m3.mirdb.v2023.2.Mm.symbols.gmt' : ''
M5       = M5       == 'true' && species == 'mouse' ? ' m5.all.v2023.2.Mm.symbols.gmt'   : ''    
M5_GO    = M5_GO    == 'true' && species == 'mouse' ? ' m5.go.v2023.2.Mm.symbols.gmt'    : '' 
M5_GO_BP = M5_GO_BP == 'true' && species == 'mouse' ? ' m5.go.bp.v2023.2.Mm.symbols.gmt' : ''
M5_GO_CC = M5_GO_CC == 'true' && species == 'mouse' ? ' m5.go.cc.v2023.2.Mm.symbols.gmt' : ''
M5_GO_MF = M5_GO_MF == 'true' && species == 'mouse' ? ' m5.go.mf.v2023.2.Mm.symbols.gmt' : ''
M5_MPT   = M5_MPT   == 'true' && species == 'mouse' ? ' m5.mpt.v2023.2.Mm.symbols.gmt'   : ''
M8       = M8       == 'true' && species == 'mouse' ? ' m8.all.v2023.2.Mm.symbols.gmt'   : ''

gmt_list = H + C1 + C2 + C2_CGP + C2_CP + C3_MIR + C3_TFT + C4 + C4_CGN + C4_CM + C5 + C5_GO + C5_GO_BP + C5_GO_CC + C5_GO_MF + C5_HPO + C6 + C7 + C8 + MH + M1 + M2 + M2_CGP + M2_CP + M3_GTRD + M3_miRDB + M5 + M5_GO + M5_GO_BP + M5_GO_CC + M5_GO_MF + M5_MPT + M8

if (gmt_list.equals("")){
	gmt_list = 'h.all.v2023.2.Hs.symbols.gmt'
}

"""
prepare_GSEA.py \
--input ${input} --species ${species} --event-column ${event_column} --fold-change-column ${fold_change_column} \
--GMT-key /data/gmt_key.txt --GMT-source /data --GMT-list ${gmt_list} --minSize ${minSize} --maxSize ${maxSize} \
--NES ${nes_significance_cutoff} --pvalue ${padj_significance_cutoff} \
--seed ${seed} --threads ${task.cpus} --postfix '_gsea${postfix}'

cp -R outputs GSEA/

mkdir GSEA_reports
cp -R GSEA GSEA_reports/
"""

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 4 
}
//* platform
//* platform
//* autofill

process DE_module_LimmaVoom_Analysis {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /(.*.Rmd|.*.html|inputs|outputs|GSEA)$/) "limmaVoom/$filename"}
input:
 path DE_reports
 path GSEA_reports

output:
 path "{*.Rmd,*.html,inputs,outputs,GSEA}"  ,emit:g37_39_outputDir00 

container 'quay.io/viascientific/de_module:2.1.0'

script:

"""
mv DE_reports/* .

if [ -d "GSEA_reports" ]; then
    mv GSEA_reports/* .
fi
"""
}


workflow {


if (!((params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal)){
g_15_0_g_31.set{g_31_reads01_g14_31}
g_31_log_file11 = Channel.empty()
} else {

Adapter_Removal(g_15_0_g_31,g_16_1_g_31)
g_31_reads01_g14_31 = Adapter_Removal.out.g_31_reads01_g14_31
g_31_log_file11 = Adapter_Removal.out.g_31_log_file11
}


FastQC(g_15_0_g_33,g_16_1_g_33)
g_33_FastQCout00 = FastQC.out.g_33_FastQCout00


featureCounts_Prep(g_36_0_g_35)
g_35_run_parameters02_g_27 = featureCounts_Prep.out.g_35_run_parameters02_g_27


STAR_Module_Check_Build_STAR_Index(g_30_0_g14_21,g_12_1_g14_21)
g14_21_starIndex00_g14_26 = STAR_Module_Check_Build_STAR_Index.out.g14_21_starIndex00_g14_26

g14_21_starIndex00_g14_26= g14_21_starIndex00_g14_26.ifEmpty(ch_empty_file_1) 


if (!((params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR)){
g14_21_starIndex00_g14_26.set{g14_26_starIndex02_g14_31}
} else {

STAR_Module_check_STAR_files(g14_21_starIndex00_g14_26)
g14_26_starIndex02_g14_31 = STAR_Module_check_STAR_files.out.g14_26_starIndex02_g14_31
}


STAR_Module_Map_STAR(g_16_0_g14_31,g_31_reads01_g14_31,g14_26_starIndex02_g14_31)
g14_31_outputFileOut00_g14_18 = STAR_Module_Map_STAR.out.g14_31_outputFileOut00_g14_18
g14_31_outputFileTxt11 = STAR_Module_Map_STAR.out.g14_31_outputFileTxt11
g14_31_logOut22 = STAR_Module_Map_STAR.out.g14_31_logOut22
g14_31_mapped_reads30_g14_30 = STAR_Module_Map_STAR.out.g14_31_mapped_reads30_g14_30
g14_31_outputFileTab44 = STAR_Module_Map_STAR.out.g14_31_outputFileTab44
g14_31_progressOut55 = STAR_Module_Map_STAR.out.g14_31_progressOut55
g14_31_transcriptome_bam60_g14_15 = STAR_Module_Map_STAR.out.g14_31_transcriptome_bam60_g14_15


STAR_Module_STAR_Summary(g14_31_outputFileOut00_g14_18.groupTuple())
g14_18_outputFile00_g14_11 = STAR_Module_STAR_Summary.out.g14_18_outputFile00_g14_11
g14_18_name11_g14_11 = STAR_Module_STAR_Summary.out.g14_18_name11_g14_11


STAR_Module_merge_tsv_files_with_same_header(g14_18_outputFile00_g14_11.collect(),g14_18_name11_g14_11.collect())
g14_11_outputFileTSV00 = STAR_Module_merge_tsv_files_with_same_header.out.g14_11_outputFileTSV00


STAR_Module_Merge_Bam_and_create_sense_antisense(g14_31_mapped_reads30_g14_30.groupTuple(),g_16_1_g14_30)
g14_30_bam_index00 = STAR_Module_Merge_Bam_and_create_sense_antisense.out.g14_30_bam_index00
g14_30_bamFile10_g_27 = STAR_Module_Merge_Bam_and_create_sense_antisense.out.g14_30_bamFile10_g_27


featureCounts(g14_30_bamFile10_g_27,g_16_1_g_27,g_35_run_parameters02_g_27,g_12_3_g_27)
g_27_outputFileTSV00_g_28 = featureCounts.out.g_27_outputFileTSV00_g_28


featureCounts_summary(g_27_outputFileTSV00_g_28.collect())
g_28_outputFile00_g37_24 = featureCounts_summary.out.g_28_outputFile00_g37_24
(g_28_outputFile00_g37_42) = [g_28_outputFile00_g37_24]
g_28_outFileTSV11 = featureCounts_summary.out.g_28_outFileTSV11


STAR_Module_merge_transcriptome_bam(g14_31_transcriptome_bam60_g14_15.groupTuple())
g14_15_merged_bams00 = STAR_Module_merge_transcriptome_bam.out.g14_15_merged_bams00
g14_15_bam_index11 = STAR_Module_merge_transcriptome_bam.out.g14_15_bam_index11
g14_15_sorted_bam22 = STAR_Module_merge_transcriptome_bam.out.g14_15_sorted_bam22



DE_module_Prepare_DESeq2(g_28_outputFile00_g37_24,g_3_1_g37_24,g_4_2_g37_24,g_21_3_g37_24)
g37_24_outputFile00_g37_37 = DE_module_Prepare_DESeq2.out.g37_24_outputFile00_g37_37
g37_24_postfix10_g37_33 = DE_module_Prepare_DESeq2.out.g37_24_postfix10_g37_33
g37_24_outputFile21_g37_33 = DE_module_Prepare_DESeq2.out.g37_24_outputFile21_g37_33

g37_24_postfix10_g37_33= g37_24_postfix10_g37_33.ifEmpty("") 


DE_module_Prepare_GSEA_DESeq2(g37_24_postfix10_g37_33,g37_24_outputFile21_g37_33,g_39_2_g37_33)
g37_33_outputFile01_g37_37 = DE_module_Prepare_GSEA_DESeq2.out.g37_33_outputFile01_g37_37
g37_33_outputFile11 = DE_module_Prepare_GSEA_DESeq2.out.g37_33_outputFile11

g37_33_outputFile01_g37_37= g37_33_outputFile01_g37_37.ifEmpty(ch_empty_file_1) 


DE_module_DESeq2_Analysis(g37_24_outputFile00_g37_37,g37_33_outputFile01_g37_37)
g37_37_outputDir00 = DE_module_DESeq2_Analysis.out.g37_37_outputDir00



DE_module_Prepare_LimmaVoom(g_28_outputFile00_g37_42,g_3_1_g37_42,g_4_2_g37_42,g_38_3_g37_42)
g37_42_outputFile00_g37_39 = DE_module_Prepare_LimmaVoom.out.g37_42_outputFile00_g37_39
g37_42_postfix10_g37_41 = DE_module_Prepare_LimmaVoom.out.g37_42_postfix10_g37_41
g37_42_outputFile21_g37_41 = DE_module_Prepare_LimmaVoom.out.g37_42_outputFile21_g37_41

g37_42_postfix10_g37_41= g37_42_postfix10_g37_41.ifEmpty("") 


DE_module_Prepare_GSEA_LimmaVoom(g37_42_postfix10_g37_41,g37_42_outputFile21_g37_41,g_40_2_g37_41)
g37_41_outputFile01_g37_39 = DE_module_Prepare_GSEA_LimmaVoom.out.g37_41_outputFile01_g37_39
g37_41_outputFile11 = DE_module_Prepare_GSEA_LimmaVoom.out.g37_41_outputFile11

g37_41_outputFile01_g37_39= g37_41_outputFile01_g37_39.ifEmpty(ch_empty_file_1) 


DE_module_LimmaVoom_Analysis(g37_42_outputFile00_g37_39,g37_41_outputFile01_g37_39)
g37_39_outputDir00 = DE_module_LimmaVoom_Analysis.out.g37_39_outputDir00


}

workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
