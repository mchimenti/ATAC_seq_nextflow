#!/usr/bin/env nextflow

genome = params.genome
bwa = params.bwa

/*
 * this is a queue channel, outputs a tuple with sample name, file names
 */

Channel
	.fromFilePairs( params.reads, flat: true )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.into { reads_ch1; reads_ch2 }

process fastqc {
	publishDir "qc_report"
	tag "sample: ${pair_id}"

	input:
	set val(pair_id), file(fastq_R1), file(fastq_R2) from reads_ch1
	
	output:
	file(output_R1_html)
        file(output_R1_zip)
        file(output_R2_html)
        file(output_R2_zip) into fastqc_ch
	

	script:
	output_R1_html = "${fastq_R1}".replaceFirst(/.fastq$/, "_fastqc.html")
        output_R1_zip = "${fastq_R1}".replaceFirst(/.fastq$/, "_fastqc.zip")
        output_R2_html = "${fastq_R2}".replaceFirst(/.fastq$/, "_fastqc.html")
        output_R2_zip = "${fastq_R2}".replaceFirst(/.fastq$/, "_fastqc.zip")
        """
        fastqc -o . "${fastq_R1}"
        fastqc -o . "${fastq_R2}"
        """
}

process ngmerge {

	maxForks 8
	tag "sample: $pair_id"

	input: 
	set val(pair_id), file(read1), file(read2) from reads_ch2
	
	output: 
	set pair_id, file("${pair_id}_trim_1.fastq"), file("${pair_id}_trim_2.fastq") into trim_read_ch
	
	script:
	"""
	~/NGmerge/NGmerge -a -v -y -u 41 -1 ${read1} -2 ${read2} -o ${pair_id}_trim 
	"""

}

process bwa {
	publishDir "aligned_BAMs"
	maxForks 8
	tag "Sample: $pair_id"

	input: 
	set val(pair_id), file(read1), file(read2) from trim_read_ch

	output:
	set pair_id, file("${pair_id}_sort.bam") into bam_ch	

	script:
	"""
	$bwa mem -c 250 -t 1 -v 3 ${genome} ${read1} ${read2} | samtools view -u - | samtools sort -n -o ${pair_id}_sort.bam  
	""" 

}

process genrich {

	publishDir "genrich_peaks"


	input:
	set val(pair_id), file("${pair_id}_sort.bam") from bam_ch

	output:
	file("${pair_id}_genrich.narrowPeak")

	script:
	"""
	~/Genrich/Genrich -t ${pair_id}_sort.bam -o ${pair_id}_genrich.narrowPeak -v -j -y -r -e MT,GL000191.1,GL000192.1,GL000193.1,GL000195.1,GL000199.1,GL000205.1,GL000206.1,GL000208.1,GL000212.1,GL000214.1,GL000216.1,GL000217.1,GL000219.1,GL000220.1,GL000222.1,GL000223.1,GL000224.1,GL000225.1,GL000226.1,GL000228.1,GL000235.1,GL000243.1
	"""

}


