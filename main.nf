$HOSTNAME = ""
params.outdir = 'results'  

evaluate(new File("${params.projectDir}/nextflow_header.config"))
//params.metadata.metadata = "${params.projectDir}/tools.json"
if (!params.reads){params.reads = ""} 
if (!params.mate_p){params.mate_p = ""} 
if (!params.mate_s){params.mate_s = ""} 

if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_0_reads_g_2}
 } else {  
	g_0_reads_g_2 = Channel.empty()
 }

Channel.value(params.mate_p).into{g_1_mate_g_2;g_1_mate_g11_15;g_1_mate_g11_19;g_1_mate_g11_12}
Channel.value(params.mate_s).into{g_23_mate_g3_0;g_23_mate_g3_5;g_23_mate_g3_7;g_23_mate_g15_9;g_23_mate_g15_12;g_23_mate_g15_11}


process unizp {

input:
 set val(name),file(reads) from g_0_reads_g_2
 val mate from g_1_mate_g_2

output:
 set val(name),file("*.fastq")  into g_2_reads0_g11_12

script:

if(mate=="pair"){
	readArray = reads.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]
	
	"""
	case "$R1" in
	*.gz | *.tgz ) 
	        gunzip -c $R1 > R1.fastq
	        ;;
	*)
	        cp $R1 ./R1.fastq
	        echo "$R1 not gzipped"
	        ;;
	esac
	
	case "$R2" in
	*.gz | *.tgz ) 
	        gunzip -c $R2 > R2.fastq
	        ;;
	*)
	        cp $R2 ./R2.fastq
	        echo "$R2 not gzipped"
	        ;;
	esac
	"""
}else{
	"""
	case "$reads" in
	*.gz | *.tgz ) 
	        gunzip -c $reads > R1.fastq
	        ;;
	*)
	        cp $reads ./R1.fastq
	        echo "$reads not gzipped"
	        ;;
	esac
	"""
}
}


process Assemble_pairs_assemble_pairs {

input:
 set val(name),file(reads) from g_2_reads0_g11_12
 val mate from g_1_mate_g11_12

output:
 set val(name),file("*_assemble-pass.f*")  into g11_12_reads0_g3_0
 set val(name),file("AP_*")  into g11_12_logFile1_g11_15
 set val(name),file("*_assemble-fail.f*") optional true  into g11_12_reads_failed22
 set val(name),file("out*")  into g11_12_logFile33

script:
method = params.Assemble_pairs_assemble_pairs.method
coord = params.Assemble_pairs_assemble_pairs.coord
rc = params.Assemble_pairs_assemble_pairs.rc
head_fields_R1 = params.Assemble_pairs_assemble_pairs.head_fields_R1
head_fields_R2 = params.Assemble_pairs_assemble_pairs.head_fields_R2
failed = params.Assemble_pairs_assemble_pairs.failed
fasta = params.Assemble_pairs_assemble_pairs.fasta
nproc = params.Assemble_pairs_assemble_pairs.nproc
alpha = params.Assemble_pairs_assemble_pairs.alpha
maxerror = params.Assemble_pairs_assemble_pairs.maxerror
minlen = params.Assemble_pairs_assemble_pairs.minlen
maxlen = params.Assemble_pairs_assemble_pairs.maxlen
scanrev = params.Assemble_pairs_assemble_pairs.scanrev
minident = params.Assemble_pairs_assemble_pairs.minident
evalue = params.Assemble_pairs_assemble_pairs.evalue
maxhits = params.Assemble_pairs_assemble_pairs.maxhits
fill = params.Assemble_pairs_assemble_pairs.fill
aligner = params.Assemble_pairs_assemble_pairs.aligner
// align_exec = params.Assemble_pairs_assemble_pairs.// align_exec
// dbexec = params.Assemble_pairs_assemble_pairs.// dbexec
gap = params.Assemble_pairs_assemble_pairs.gap
usearch_version = params.Assemble_pairs_assemble_pairs.usearch_version
assemble_reference = params.Assemble_pairs_assemble_pairs.assemble_reference
head_seqeunce_file = params.Assemble_pairs_assemble_pairs.head_seqeunce_file
//* @style @condition:{method="align",alpha,maxerror,minlen,maxlen,scanrev}, {method="sequential",alpha,maxerror,minlen,maxlen,scanrev,ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="reference",ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="join",gap} @multicolumn:{method,coord,rc,head_fields_R1,head_fields_R2,failed,nrpoc,usearch_version},{alpha,maxerror,minlen,maxlen,scanrev}, {ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec}, {gap} 

// args
coord = "--coord ${coord}"
rc = "--rc ${rc}"
head_fields_R1 = (head_fields_R1!="") ? "--1f ${head_fields_R1}" : ""
head_fields_R2 = (head_fields_R2!="") ? "--2f ${head_fields_R2}" : ""
failed = (failed=="false") ? "" : "--failed"
fasta = (fasta=="false") ? "" : "--fasta"
nproc = "--nproc ${nproc}"

scanrev = (scanrev=="false") ? "" : "--scanrev"
fill = (fill=="false") ? "" : "--fill"

// align_exec = (align_exec!="") ? "--exec ${align_exec}" : ""
// dbexec = (dbexec!="") ? "--dbexec ${dbexec}" : ""


ref_file = (assemble_reference!='') ? "-r ${assemble_reference}" : ""



args = ""

if(method=="align"){
	args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev}"
}else{
	if(method=="sequential"){
		args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev} ${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
	}else{
		if(method=="reference"){
			args = "${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
		}else{
			args = "--gap ${gap}"
		}
	}
}


readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	
	if(R1.contains(""+head_seqeunce_file)){
		R1 = readArray[0]
		R2 = readArray[1]
	}else{
		R2 = readArray[0]
		R1 = readArray[1]
	}
	
	"""
	if [ "${method}" != "align" ]; then
		if  [ "${aligner}" == "usearch" ]; then
			wget -q --show-progress --no-check-certificate https://drive5.com/downloads/usearch${usearch_version}_i86linux32.gz
			gunzip usearch${usearch_version}_i86linux32.gz
			chmod +x usearch${usearch_version}_i86linux32
			mv usearch${usearch_version}_i86linux32 /usr/local/bin/usearch2
			align_exec="--exec /usr/local/bin/usearch2"
			dbexec="--dbexec /usr/local/bin/usearch2"
		else
			align_exec="--exec /usr/local/bin/blastn"
			dbexec="--dbexec /usr/local/bin/makeblastdb"
		fi
	else
		align_exec=""
		dbexec=""
	fi

	AssemblePairs.py ${method} -1 ${R1} -2 ${R2} ${coord} ${rc} ${head_fields_R1} ${head_fields_R2} ${args} \$align_exec \$dbexec ${fasta} ${failed} --log AP_${name}.log ${nproc}  2>&1 | tee out_${R1}_AP.log
	"""

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}

}


process Filter_Sequence_Quality_filter_seq_quality {

input:
 set val(name),file(reads) from g11_12_reads0_g3_0
 val mate from g_23_mate_g3_0

output:
 set val(name), file("*_${method}-pass.fast*")  into g3_0_reads0_g15_11
 set val(name), file("FS_*")  into g3_0_logFile1_g3_5
 set val(name), file("*_${method}-fail.fast*") optional true  into g3_0_reads22
 set val(name),file("out*") optional true  into g3_0_logFile33

script:
method = params.Filter_Sequence_Quality_filter_seq_quality.method
nproc = params.Filter_Sequence_Quality_filter_seq_quality.nproc
q = params.Filter_Sequence_Quality_filter_seq_quality.q
n_length = params.Filter_Sequence_Quality_filter_seq_quality.n_length
n_missing = params.Filter_Sequence_Quality_filter_seq_quality.n_missing
fasta = params.Filter_Sequence_Quality_filter_seq_quality.fasta
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="missing"){
	q = ""
	n_length = ""
	n_missing = "-n ${n_missing}"
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
	}else{
		q = "-q ${q}"
		n_length = ""
		n_missing = ""
	}
}

readArray = reads.toString().split(' ')	

fasta = (fasta=="true") ? "--fasta" : ""

if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R1_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R2_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}


}


process Filter_Sequence_Quality_parse_log_FS {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*table.tab$/) "FSQ_log_table/$filename"}
input:
 set val(name), file(log_file) from g3_0_logFile1_g3_5
 val mate from g_23_mate_g3_5

output:
 file "*table.tab"  into g3_5_logFile0_g3_7, g3_5_logFile0_g3_16

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f ID QUALITY
"""

}


process Filter_Sequence_Quality_report_filter_Seq_Quality {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.rmd$/) "FSQ_report/$filename"}
input:
 val mate from g_23_mate_g3_7
 file log_files from g3_5_logFile0_g3_7

output:
 file "*.rmd"  into g3_7_rMarkdown0_g3_16


shell:

if(mate=="pair"){
	readArray = log_files.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]

	name = R1 - "_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read 1", "Read 2")
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("quality", 
	        paste("Mean Phred quality scores for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The dotted line indicates the average quality score under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	quality_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	quality_log_2 <- loadLogTable(file.path(".", "!{R2}"))
	```
	
	# Quality Scores
	
	Quality filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq tool remove reads with low mean Phred quality scores. 
	Phred quality scores are assigned to each nucleotide base call in automated 
	sequencer traces. The quality score (`Q`) of a base call is logarithmically 
	related to the probability that a base call is incorrect (`P`): 
	$Q = -10 log_{10} P$. For example, a base call with `Q=30` is incorrectly 
	assigned 1 in 1000 times. The most commonly used approach is to remove read 
	with average `Q` below 20.
	
	```{r, echo=FALSE}
	plotFilterSeq(quality_log_1, quality_log_2, titles=plot_titles, sizing="figure")
	```
	
	`r figures("quality")`
		
	EOF
	
	open OUT, ">FSQ_!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = log_files.toString().split(' ')
	R1 = readArray[0]
	name = R1 - "_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read")#params$quality_titles
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("quality", 
	        paste("Mean Phred quality scores for",  plot_titles[1],
	              "The dotted line indicates the average quality score under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	quality_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	```
	
	# Quality Scores
	
	Quality filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq tool remove reads with low mean Phred quality scores. 
	Phred quality scores are assigned to each nucleotide base call in automated 
	sequencer traces. The quality score (`Q`) of a base call is logarithmically 
	related to the probability that a base call is incorrect (`P`): 
	$Q = -10 log_{10} P$. For example, a base call with `Q=30` is incorrectly 
	assigned 1 in 1000 times. The most commonly used approach is to remove read 
	with average `Q` below 20.
	
	```{r, echo=FALSE}
	plotFilterSeq(quality_log_1, titles=plot_titles[1], sizing="figure")
	```
	
	`r figures("quality")`
	
	EOF
	
	open OUT, ">FSQ_!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Filter_Sequence_Quality_presto_render_rmarkdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.html$/) "FSQ_html_report/$filename"}
input:
 file rmk from g3_7_rMarkdown0_g3_16
 file log_file from g3_5_logFile0_g3_16

output:
 file "*.html" optional true  into g3_16_outputFileHTML00
 file "*csv" optional true  into g3_16_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Assemble_pairs_parse_log_AP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*table.tab$/) "AP_log_table/$filename"}
input:
 set val(name),file(log_file) from g11_12_logFile1_g11_15
 val mate from g_1_mate_g11_15

output:
 file "*table.tab"  into g11_15_logFile0_g11_25, g11_15_logFile0_g11_19

script:
field_to_parse = params.Assemble_pairs_parse_log_AP.field_to_parse
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ${field_to_parse}
"""


}


process Assemble_pairs_report_assemble_pairs {

input:
 file log_files from g11_15_logFile0_g11_19
 val matee from g_1_mate_g11_19

output:
 file "*.rmd"  into g11_19_rMarkdown0_g11_25



shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')
	assemble = readArray[0]
	name = assemble-"_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("assemble_length", "Histogram showing the distribution assembled sequence lengths in 
	                            nucleotides for the Align step (top) and Reference step (bottom).")
	figures("assemble_overlap", "Histogram showing the distribution of overlapping nucleotides between 
	                             mate-pairs for the Align step (top) and Reference step (bottom).
	                             Negative values for overlap indicate non-overlapping mate-pairs
	                             with the negative value being the number of gap characters between
	                             the ends of the two mate-pairs.")
	figures("assemble_error", "Histograms showing the distribution of paired-end assembly error 
	                           rates for the Align step (top) and identity to the reference germline 
	                           for the Reference step (bottom).")
	figures("assemble_pvalue", "Histograms showing the distribution of significance scores for 
	                            paired-end assemblies. P-values for the Align mode are shown in the top
	                            panel. E-values from the Reference step's alignment against the 
	                            germline sequences are shown in the bottom panel for both input files
	                            separately.")
	```
	
	```{r, echo=FALSE, warning=FALSE}
	assemble_log <- loadLogTable(file.path(".", "!{assemble}"))
	
	# Subset to align and reference logs
	align_fields <- c("ERROR", "PVALUE")
	ref_fields <- c("REFID", "GAP", "EVALUE1", "EVALUE2", "IDENTITY")
	align_log <- assemble_log[!is.na(assemble_log$ERROR), !(names(assemble_log) %in% ref_fields)]
	ref_log <- assemble_log[!is.na(assemble_log$REFID), !(names(assemble_log) %in% align_fields)]
	
	# Build log set
	assemble_list <- list()
	if (nrow(align_log) > 0) { assemble_list[["Align"]] <- align_log }
	if (nrow(ref_log) > 0) { assemble_list[["Reference"]] <- ref_log }
	plot_titles <- names(assemble_list)
	```
	
	# Paired-End Assembly
	
	Assembly of paired-end reads is performed using the AssemblePairs tool which 
	determines the read overlap in two steps. First, de novo assembly is attempted 
	using an exhaustive approach to identify all possible overlaps between the 
	two reads with alignment error rates and p-values below user-defined thresholds. 
	This method is denoted as the `Align` method in the following figures. 
	Second, those reads failing the first stage of de novo assembly are then 
	mapped to the V-region reference sequences to create a full length sequence, 
	padding with Ns, for any amplicons that have insufficient overlap for 
	de novo assembly. This second stage is referred to as the `Reference` step in the
	figures below.
	
	## Assembled sequence lengths
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="length", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_length")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="overlap", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_overlap")`
	
	## Alignment error rates and significance
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="error", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_error")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="pvalue", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```

	`r figures("assemble_pvalue")`

	EOF
	
	open OUT, ">AP_!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}
}


process Assemble_pairs_presto_render_rmarkdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.html$/) "AP_report/$filename"}
input:
 file rmk from g11_19_rMarkdown0_g11_25
 file log_file from g11_15_logFile0_g11_25

output:
 file "*.html" optional true  into g11_25_outputFileHTML00
 file "*csv" optional true  into g11_25_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Mask_Primer_align_MaskPrimers {

input:
 val mate from g_23_mate_g15_11
 set val(name),file(reads) from g3_0_reads0_g15_11

output:
 set val(name), file("*_primers-pass.fast*") optional true  into g15_11_reads0_g17_16
 set val(name), file("*_primers-fail.fast*") optional true  into g15_11_reads_failed11
 set val(name), file("MP_*")  into g15_11_logFile2_g15_9
 set val(name),file("out*")  into g15_11_logFile33

script:
method = params.Mask_Primer_align_MaskPrimers.method
barcode_field = params.Mask_Primer_align_MaskPrimers.barcode_field
primer_field = params.Mask_Primer_align_MaskPrimers.primer_field
barcode = params.Mask_Primer_align_MaskPrimers.barcode
revpr = params.Mask_Primer_align_MaskPrimers.revpr
mode = params.Mask_Primer_align_MaskPrimers.mode
failed = params.Mask_Primer_align_MaskPrimers.failed
fasta = params.Mask_Primer_align_MaskPrimers.fasta
nproc = params.Mask_Primer_align_MaskPrimers.nproc
maxerror = params.Mask_Primer_align_MaskPrimers.maxerror
umi_length = params.Mask_Primer_align_MaskPrimers.umi_length
start = params.Mask_Primer_align_MaskPrimers.start
extract_length = params.Mask_Primer_align_MaskPrimers.extract_length
maxlen = params.Mask_Primer_align_MaskPrimers.maxlen
skiprc = params.Mask_Primer_align_MaskPrimers.skiprc
R1_primers = params.Mask_Primer_align_MaskPrimers.R1_primers
R2_primers = params.Mask_Primer_align_MaskPrimers.R2_primers
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""
fasta = (fasta=="true") ? "--fasta" : ""
def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    PRIMER_FIELD = "${pf}"
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk} ${pf} ${bf}")
    
    
}}

readArray = reads.toString().split(' ')
if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
  


	R1 = readArray[0]
	R2 = readArray[1]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	
	"""
	
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} ${fasta} 2>&1 | tee -a out_${name}_MP.log & \
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} ${fasta} 2>&1 | tee -a out_${name}_MP.log & \
	wait
	"""
}else{
	args_1 = args_values[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	R1 = readArray[0]

	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} ${fasta} 2>&1 | tee -a out_${name}_MP.log
	"""
}

}


process collapse_sequences_collapse_seq {

input:
 set val(name), file(reads) from g15_11_reads0_g17_16

output:
 set val(name),  file("*_collapse-unique.fast*")  into g17_16_reads0_g18_20
 set val(name),  file("*_collapse-duplicate.fast*") optional true  into g17_16_reads_duplicate11
 set val(name),  file("*_collapse-undetermined.fast*") optional true  into g17_16_reads_undetermined22
 file "CS_*"  into g17_16_logFile33

script:
max_missing = params.collapse_sequences_collapse_seq.max_missing
inner = params.collapse_sequences_collapse_seq.inner
fasta = params.collapse_sequences_collapse_seq.fasta
act = params.collapse_sequences_collapse_seq.act
uf = params.collapse_sequences_collapse_seq.uf
cf = params.collapse_sequences_collapse_seq.cf
nproc = params.collapse_sequences_collapse_seq.nproc
failed = params.collapse_sequences_collapse_seq.failed

inner = (inner=="true") ? "--inner" : ""
fasta = (fasta=="true") ? "--fasta" : ""
act = (act=="none") ? "" : "--act ${act}"
cf = (cf=="") ? "" : "--cf ${cf}"
uf = (uf=="") ? "" : "--uf ${uf}"
failed = (failed=="false") ? "" : "--failed"

"""
CollapseSeq.py -s ${reads} -n ${max_missing} ${fasta} ${inner} ${uf} ${cf} ${act} --log CS_${name}.log ${failed}
"""

}


process split_sequences_split_seq {

input:
 set val(name),file(reads) from g17_16_reads0_g18_20

output:
 set val(name), file("*_atleast-*.fast*")  into g18_20_reads0_g_27, g18_20_reads0_g19_15

script:
field = params.split_sequences_split_seq.field
num = params.split_sequences_split_seq.num
fasta = params.split_sequences_split_seq.fasta

readArray = reads.toString()

if(num!=0){
	num = " --num ${num}"
}else{
	num = ""
}

fasta = (fasta=="false") ? "" : "--fasta"

"""
SplitSeq.py group -s ${readArray} -f ${field} ${num} ${fasta}
"""

}


process Split_TCR_chains {

input:
 set val(name), file(reads) from g18_20_reads0_g_27

output:
 set val(name),  file("TRA.fastq")  into g_27_reads0_g_28
 set val(name),  file("TRB.fastq")  into g_27_reads1_g_29
 set val(name),  file("UNKNOWN.fastq")  into g_27_reads22

script:
nproc = params.Split_TCR_chains.nproc

readArray = reads.toString().split(' ')
R1 = readArray[0]

// Define output filenames
def out_TRA = "TRA.fastq"
def out_TRB = "TRB.fastq"
def out_UNKNOWN = "UNKNOWN.fastq"

"""
${fasta} ${R1} | awk '
    BEGIN {
        OFS="\\n";
    }
    {
        if (NR % 4 == 1) {
            header = \$0;
            if (header ~ /PRIMER=TRA/) {
                file = "${out_TRA}";
            }
            else if (header ~ /PRIMER=TRB/) {
                file = "${out_TRB}";
            }
            else {
                file = "${out_UNKNOWN}";
            }
            print header > file;
            getline seq;
            print seq > file;
            getline plus;
            print plus > file;
            getline qual;
            print qual > file;
        }
    }
'
"""


}


process vdjbase_input_trb {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${chain}$/) "reads/$filename"}
input:
 set val(name),file(reads) from g_27_reads1_g_29

output:
 file "${chain}"  into g_29_germlineDb00

script:
chain = params.vdjbase_input_trb.chain

"""
mkdir ${chain}
mv ${reads} ${chain}/${name}.fasta
"""

}


process vdjbase_input_tra {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${chain}$/) "reads/$filename"}
input:
 set val(name),file(reads) from g_27_reads0_g_28

output:
 file "${chain}"  into g_28_germlineDb00

script:
chain = params.vdjbase_input_tra.chain

"""
mkdir ${chain}
mv ${reads} ${chain}/${name}.fasta
"""

}


process Parse_header_table_parse_headers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*${out}$/) "parse_header_table/$filename"}
input:
 set val(name), file(reads) from g18_20_reads0_g19_15

output:
 set val(name),file("*${out}")  into g19_15_reads00

script:
method = params.Parse_header_table_parse_headers.method
act = params.Parse_header_table_parse_headers.act
args = params.Parse_header_table_parse_headers.args


if(method=="collapse" || method=="copy" || method=="rename" || method=="merge"){
	out="_reheader.fastq"
	act = (act=="none") ? "" : "--act ${act}"
	"""
	ParseHeaders.py  ${method} -s ${reads} ${args} ${act}
	"""
}else{
	if(method=="table"){
			out=".tab"
			"""
			ParseHeaders.py ${method} -s ${reads} ${args}
			"""	
	}else{
		out="_reheader.fastq"
		"""
		ParseHeaders.py ${method} -s ${reads} ${args}
		"""		
	}
}


}


process Mask_Primer_align_parse_log_MP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*table.tab$/) "MP_align_log/$filename"}
input:
 val mate from g_23_mate_g15_9
 set val(name), file(log_file) from g15_11_logFile2_g15_9

output:
 file "*table.tab"  into g15_9_logFile0_g15_12, g15_9_logFile0_g15_19

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process Mask_Primer_align_try_report_maskprimer {

input:
 file primers from g15_9_logFile0_g15_12
 val mate from g_23_mate_g15_12

output:
 file "*.rmd"  into g15_12_rMarkdown0_g15_19


shell:

if(mate=="pair"){
	readArray = primers.toString().split(' ')	
	primers_1 = readArray[0]
	primers_2 = readArray[1]
	name = primers_1 - "_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read 1", "Read 2")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom),",
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers_1}"))
	primer_log_2 <- loadLogTable(file.path(".", "!{primers_2}"))
	
	primer_log1_error <- any(is.na(primer_log_1[['ERROR']]))
	primer_log2_error<- any(is.na(primer_log_2[['ERROR']]))
	
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	if(!primer_log1_error && !primer_log2_error)
		plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	if(!primer_log1_error && !primer_log2_error)
		plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	# check the error column exists 
	if(!primer_log1_error && !primer_log2_error)
		plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = primers.toString().split(' ')
	primers = readArray[0]
	name = primers - "_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1],
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Mask_Primer_align_presto_render_rmarkdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.html$/) "MP_align_report/$filename"}
input:
 file rmk from g15_12_rMarkdown0_g15_19
 file log_file from g15_9_logFile0_g15_19

output:
 file "*.html" optional true  into g15_19_outputFileHTML00
 file "*csv" optional true  into g15_19_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
