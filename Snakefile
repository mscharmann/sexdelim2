configfile: "config.yaml"

popmapfile = config["popmapfile"]
windowsize = int( config["windowsize"] )
genomefile = config["genomefile"]
samples_units_fqs_map = config["samples_units_fqs_map"]
regions_for_plot_bed = config["regions_for_plot_bed"]


import pandas as pd

samples_units_fqs = pd.read_table(samples_units_fqs_map, dtype=str).set_index(
	["sample", "unit", "fq1", "fq2"], drop=False)

SAMPLES = list( set(samples_units_fqs["sample"]) )


def get_sample_bams(wildcards):
	"""Get all aligned reads of given sample."""
	return expand(
		"mapped_reads_per_unit/{sample}-{unit}.sorted.bam",
		sample=wildcards.sample,
		unit=samples_units_fqs.loc[wildcards.sample].unit,
	)


def get_fastq_sample_unit(wildcards):
	"""Get fastq files of given sample-unit."""
	fastqs = samples_units_fqs.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]]
	if fastqs.fq2.isnull().values.any():
		return [ fastqs.fq1.item() ]
	return [ fastqs.fq1.item(), fastqs.fq2.item() ]


def get_fastq_sample_ALL(wildcards):
	"""Get list of fastq files of given sample including ALL units."""
	fastqs_pd = samples_units_fqs.loc[(wildcards.sample), ["fq1", "fq2"]]
	fastqs = set( fastqs_pd.fq1.tolist() + fastqs_pd.fq2.tolist() )
	fastqs_clean = list( {x for x in fastqs if pd.notna(x)} )
	return fastqs_clean



def read_popmap(popmapfile):

	SAMPLES = []
	with open(popmapfile, "r") as I:
		for line in I:
			if len(line) > 1:
				fields = line.strip("\n").split("\t")
				SAMPLES.append(fields[0])
	return(SAMPLES)


def parse_samples_reads_map (infile):

	samples_files_dict = {}
	with open(infile, "r") as I:
		I.readline() # this file has a header!
		for line in I:
			if len(line) > 1:
				fields = line.strip("\n").split("\t")
				s = fields[0]
				fq1 = fields[1]
				try:
					fq2 = fields[2]
					if len(fq2.strip()) == 0:
						fq2 = None
				except IndexError:
					fq2 = None
				try:
					samples_files_dict[s].append( (fq1, fq2) )
				except KeyError:
					samples_files_dict[s] = [ (fq1, fq2) ]

	return(samples_files_dict)


def get_fastq(wildcards):
	"""Get fastq files of given sample -- DANGER: currently at most one pair of read files per sample is handled!"""
	fqs = samples_files_dict[wildcards.sample][0]
	if not fqs[1]:
		return [fqs[0]]
	return [fqs[0], fqs[1]]



SAMPLES = read_popmap(popmapfile)
print (SAMPLES)


rule all:
	input:
		"results_raw/all.post_filter.vcf.gz",
		"results_raw/mapping_statistics_report.txt",
		"results_processed/windows.multicov.F.txt",
		"results_processed/windows.multicov.M.txt",
		"results_processed/M_F_norm_coverage_ratio_log2.bed.txt",
		"results_processed/plink.assoc_results.significant.window.bed.txt",
		"results_raw/plink.assoc_results.assoc.raw.txt.gz",
		"results_processed/LD_r2_averaged_per_window.bed.txt",
		"results_processed/kmers.F_specific.per_window.bed.txt",
		"results_processed/kmers.M_specific.per_window.bed.txt",
		"results_processed/gametolog_candidate_alleles.ZW_divergence.windows.bed.txt",
		"results_processed/gametolog_candidate_alleles.XY_divergence.windows.bed.txt",
		"results_processed/MvF.dxy.bed.txt",
		"results_processed/MvF.fst.bed.txt",
		"results_processed/F.pi.bed.txt",
		"results_processed/M.pi.bed.txt",
		"results_processed/MvF.netdiv_M.bed.txt",
		"results_processed/MvF.netdiv_F.bed.txt",
		"results_processed/het_males.bed.txt",
		"results_processed/het_females.bed.txt",
		"results_processed/allstats.txt",
		"results_processed/region.stats.txt",
		"results_processed/region.stats.plots.pdf"



rule bwa_idx:
	input:
	   genomefile
	output:
		"bwa_index/reference.fa.bwt",
		"bwa_index/reference.fa"
	shell:
		"""
		if [[ ! $( grep ">" {input} ) =~ "|" ]]; then
			mkdir -p bwa_index
			cp {input} bwa_index/reference.fa
			cd bwa_index
			bwa index reference.fa
		else
			echo "refusing to run, fasta headers contain pipe '|' character, dying"
		fi
		"""

rule bwa_map:
	input:
		fa="bwa_index/reference.fa",
		gidx="bwa_index/reference.fa.bwt",
		reads=get_fastq_sample_unit
	output:
		temp("mapped_reads_per_unit/{sample}-{unit}.bam")
	threads: 24
	run:
		if len(input.reads) == 2: # paired-end!
			shell("""
				# filtering alignments to be primary (i.e. each read only once; if multiple locations equally possible than a random one is chosen):
				# -F 256 == -F 0x0100 == NOT not primary alignment
				# filtering alignments to be NOT supplementary (supplementary: sections of read map to discontinuous coordinates, e.g. across an inversion breakpoint..):
				# -F 2048 == -F 0x800 == NOT supplementary alignment
				# sum of the bit flags: 2304 => filters against BOTH non-primary and supplementary alignments; verified with samtools flagstat
				# filtering against multi-mapping alignments (which have MAPQ=0): -q 1
				# OK, now actually retain the supplementary ones,
				# NOT reduce seed length, mismatch and clipping penalities, lower min score to output in an attemt to map more erroneous / divergent reads: -k 18 -B 3 -L 4 -T 20
				bwa mem -t {threads} -a {input.fa} {input.reads[0]} {input.reads[1]} -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:Illumina" | samtools view -F 256 -q 1 -b -@ 2 - > {output}
				""")
		else: # single-end
			shell("""
				# filtering alignments to be primary (i.e. each read only once; if multiple locations equally possible than a random one is chosen):
				# -F 256 == -F 0x0100 == NOT not primary alignment
				# filtering alignments to be NOT supplementary (supplementary: sections of read map to discontinuous coordinates, e.g. across an inversion breakpoint..):
				# -F 2048 == -F 0x800 == NOT supplementary alignment
				# -F 4 read unmapped (0x4)
				# sum of the bit flags: 2308 => filters against non-primary and supplementary alignments and unmapped
				# filtering against multi-mapping alignments (which have MAPQ=0): -q 1
				# NOT reduce seed length, mismatch and clipping penalities, lower min score to output in an attemt to map more erroneous / divergent reads: -k 18 -B 3 -L 4 -T 20
				bwa mem -t {threads} -a {input.fa} {input.reads[0]} -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:Illumina" | samtools view -F 256 -q 1 -b -@ 2 - > {output}
				""")


rule samtools_sort:
	input:
		"mapped_reads_per_unit/{sample}-{unit}.bam"
	output:
		temp( "mapped_reads_per_unit/{sample}-{unit}.sorted.bam" )
	shell:
		"""
		samtools sort -T mapped_reads_per_unit/{wildcards.sample}.{wildcards.unit} -O bam {input} > {output}
		"""

rule merge_bams_per_sample:
	input:
		bams=get_sample_bams
	output:
		"mapped_reads/{sample}.sorted.bam"
	threads: 4
	shell:
		"""
		samtools merge --threads {threads} {output} {input.bams}
		"""


rule samtools_index:
	input:
		"mapped_reads/{sample}.sorted.bam"
	output:
		"mapped_reads/{sample}.sorted.bam.bai"
	shell:
		"samtools index {input}"


# A checkpoint that shall trigger re-evaluation of the DAG
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution
# we need this because the number of output files (i.e. chunks / region files for variant calling) is not known before execution!
checkpoint split_ref_for_varcall:
	input:
		ref=genomefile,
		indexes=expand("mapped_reads/{sample}.sorted.bam.bai", sample=SAMPLES),
		samples=expand("mapped_reads/{sample}.sorted.bam", sample=SAMPLES)
	output:
		directory("varcall_chunks")
	params:
		chunksize_Mb = config["varcall_chunksize"]
	shell:
		"""
		samtools faidx {input.ref}

		echo {input.samples} | sed 's/ /\\n/g' > bamlistforsplit
		# wget https://gist.githubusercontent.com/travc/0c53df2c8eca81c3ebc36616869930ec/raw/eff3032ca7c955ca33bffd8758092e4006949c75/split_ref_by_bai_datasize.py
		python scripts/split_ref_by_bai_datasize.py -r {input.ref}.fai -L bamlistforsplit --target-data-size {params.chunksize_Mb} > varcall_regionsALL.bed
		rm bamlistforsplit
		split --lines=1 varcall_regionsALL.bed varcall_regions_ --numeric-suffixes --suffix-length=6 --additional-suffix=.bed
		mkdir -p {output}
		for i in varcall_regions_* ; do mv $i {output} ; done
		"""

rule call_variants:
	input:
		ref=genomefile,
		regions="varcall_chunks/{i}.bed"
	output:
		temp( "varcall_chunk_VCFs/{i}.bed.vcf.gz" )
	shell:
		"""
		# --max-depth 250: use at most 250 reads per input BAM file, apparently these are sampled RANDOMLY!?
		# NOT --min-MQ 20: minimum mapping quality of an alignment, otherwise skip
		# NOT --no-BAQ : do NOT re-calculate mapping quality (which involves re-aligning). Instead, will use MAPQ as stored in BAM file.
		# NOT --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [1]
		bcftools mpileup -Ou -f {input.ref} -R {input.regions} --bam-list <( ls mapped_reads/*.bam ) --max-depth 250 --min-MQ 0 --min-BQ 1 -a DP -a INFO/AD -a FORMAT/AD | bcftools call -m --skip-variants indels -Ov | bgzip -c > {output}
		"""


def merge_vcfs_input(wildcards):
	checkpoint_output = checkpoints.split_ref_for_varcall.get(**wildcards).output[0]
	# print(checkpoint_output)
	vcfs_before_filter = sorted( expand("varcall_chunk_VCFs/{i}.bed.vcf.gz", i=glob_wildcards(os.path.join(checkpoint_output, "{i}.bed")).i) )
	# print (thef)
	vcfs_after_filter = sorted( expand("varcall_chunk_VCFs_filtered/{i}.bed.vcf.gz", i=glob_wildcards(os.path.join(checkpoint_output, "{i}.bed")).i) )
	return vcfs_after_filter


def VCF_get_mean_depth_per_site_input(wildcards):
	checkpoint_output = checkpoints.split_ref_for_varcall.get(**wildcards).output[0]
	# print(checkpoint_output)
	vcfs_before_filter = sorted( expand("varcall_chunk_VCFs/{i}.bed.vcf.gz", i=glob_wildcards(os.path.join(checkpoint_output, "{i}.bed")).i) )
	return vcfs_before_filter


rule merge_filtered_vcfs:
	input:
		region_vcfs=merge_vcfs_input # refers to the function above which evaluates the checkpoint
	output:
		"results_raw/all.post_filter.vcf.gz"
	threads: 3
	shell:
		"""
		# MUST NOT USE bcftools concat: it cannot resolve POS that are non-monotonically icreasing (which ca happen at the interval boundaries)
		# the code below is only slightly modified from freebayes-parallel script: https://github.com/freebayes/freebayes/blob/master/scripts/freebayes-parallel
		# zcat input.region_vcfs | python $(which vcffirstheader) | vcfstreamsort -w 10000 | vcfuniq | bgzip -c > {output}
		# zcat alone may complain about too many arguments, so better use find -exec :
		find varcall_chunk_VCFs_filtered/*.bed.vcf.gz -type f -exec zcat {{}} \\; | python $(which vcffirstheader) | vcfstreamsort -w 10000 | vcfuniq | bgzip -c > {output}

		rm -r varcall_chunks varcall_chunk_VCFs varcall_chunk_VCFs_filtered varcall_regionsALL.bed

		"""


rule VCF_get_mean_depth_per_site:
	input:
		region_vcfs=VCF_get_mean_depth_per_site_input # function which evaluates the same checkpoint as merge_filtered_vcfs
	output:
		"results_raw/meandepthpersite.txt"
	threads: 16
	shell:
		"""
		vcf_files=$(find varcall_chunk_VCFs/*.bed.vcf.gz -type f)

		parallel -j {threads} 'zcat {{}} | grep -v "#" | cut -f1,2 > {{}}.sites' ::: ${{vcf_files}}

		find varcall_chunk_VCFs/*.bed.vcf.gz.sites -type f | xargs cat > all_pos_in_vcf.txt
		rm varcall_chunk_VCFs/*.sites

		# get a subset of 10k to 100k of the raw VCF records; iteratively keeping only 10% of lines in file
		cp all_pos_in_vcf.txt tmp1
		nlines=$( cat tmp1 | wc -l )
		echo $nlines
		while [ "$nlines " -gt 100000 ] ; do
		awk 'NR == 1 || NR % 10 == 0' tmp1 > tmp2
		mv tmp2 tmp1
		nlines=$( cat tmp1 | wc -l )
		echo $nlines
		done

		grep -v "#" tmp1 > postokeep.txt


		#parallel -j {threads} 'tabix {{}} ; tabix -h -R postokeep.txt {{}} > {{}}.subset.vcf ; rm {{}}.tbi ; vcftools --vcf {{}}.subset.vcf --site-mean-depth --stdout | cut -f3 | tail -n +2 > {{}}.mdps ; rm {{}}.subset.vcf' ::: ${{vcf_files}}

		parallel -j {threads} 'vcftools --gzvcf {{}} --positions postokeep.txt --site-mean-depth --stdout | cut -f3 | tail -n +2 > {{}}.mdps' ::: ${{vcf_files}}

		find varcall_chunk_VCFs/*.mdps -type f | xargs cat > {output}

		rm postokeep.txt all_pos_in_vcf.txt tmp1 varcall_chunk_VCFs/*.mdps

		"""


rule VCF_filter_variants_and_invariants:
	input:
		mdps="results_raw/meandepthpersite.txt",
		gzvcf="varcall_chunk_VCFs/{i}.bed.vcf.gz",
		popmap={popmapfile}
	output:
		temp( "varcall_chunk_VCFs_filtered/{i}.bed.vcf.gz" )
	params:
		MISS=config["VCF_MISS"],
		QUAL=config["VCF_QUAL"],
		MIN_DEPTH=config["VCF_MIN_DEPTH"]
	shell:
		"""
		## https://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/#comment-1469
		## => "We use ... a very vanilla filter with only depth and quality (QUAL < 20, DP < 5)"
		## => impose a minimum of QUAL, and a minimum of DEPTH
		## 	for variants, I want those two, AND a maximum DPETH cutoff
		## 	for invariants, we can only impose a minimum DEPTH and a maximum DEPTH.
		## use code from dDocent to calcualte a mean depth histogram, then find the 95th percentile: no. too complicated; there are easier ways to do this.
		## also good info:
		## 	https://speciationgenomics.github.io/filtering_vcfs/

		# prepare
		wd=DIR_{wildcards.i}
		mkdir -p varcall_chunk_VCFs_filtered/$wd
		cd varcall_chunk_VCFs_filtered/$wd

		# set filters
		# MAX_DEPTH is specific to the dataset; we take the 98th percentile of the mean depths per site:
		MAX_DEPTH=$( sort --temporary-directory=./ -n ../../{input.mdps} | awk 'BEGIN{{c=0}} length($0){{a[c]=$0;c++}}END{{p2=(c/100*2); p2=p2%1?int(p2)+1:p2; print a[c-p2-1]}}' )

		echo $MAX_DEPTH

		# variants:
		vcftools --gzvcf ../../{input.gzvcf} --mac 1 --minQ {params.QUAL} --min-meanDP {params.MIN_DEPTH} --minDP {params.MIN_DEPTH} \
		--max-meanDP $MAX_DEPTH --recode --stdout | bgzip -c > tmp.1
		tabix tmp.1

		# invariants:
		vcftools --gzvcf ../../{input.gzvcf} --max-maf 0 --min-meanDP {params.MIN_DEPTH} --minDP {params.MIN_DEPTH} \
		--max-meanDP $MAX_DEPTH --recode --stdout | bgzip -c > tmp.2
		tabix tmp.2

		# combine the two VCFs using bcftools concat:
		bcftools concat --allow-overlaps tmp.1 tmp.2 | bcftools view --exclude-uncalled --trim-alt-alleles | bgzip -c > tmp.3

		# split by M and F populations, then filter for missingness in each pop. Thus we get sites that fulfill "MISS" in at least one of either M or F populations.
		cat ../../{input.popmap} | awk '{{if($2==2) print $1}}' > fpop
		cat ../../{input.popmap} | awk '{{if($2==1) print $1}}' > mpop
		vcftools --gzvcf tmp.3 --keep fpop --max-missing {params.MISS} --recode --stdout | bgzip -c > tmp.4
		vcftools --gzvcf tmp.3 --keep mpop --max-missing {params.MISS} --recode --stdout | bgzip -c > tmp.5
		tabix tmp.4
		tabix tmp.5

		# merge M and F files again.
		bcftools merge tmp.4 tmp.5 | bgzip -c > ../../{output}

		# cleanup
		cd ../
		rm -r $wd

		"""

rule make_variant_only_VCF:
	input:
		"results_raw/all.post_filter.vcf.gz"
	output:
		temp( "results_raw/variant_sites.vcf.gz" )
	shell:
		"""
		vcftools --gzvcf {input} --mac 1 --recode --stdout | bgzip -c > {output}
		"""

rule get_coverage_in_windows:
	input:
		fa=genomefile,
		popmap={popmapfile},
		bam=expand("mapped_reads/{sample}.sorted.bam", sample=SAMPLES),
		bai=expand("mapped_reads/{sample}.sorted.bam.bai", sample=SAMPLES)
	output:
		fcov="results_processed/windows.multicov.F.txt",
		mcov="results_processed/windows.multicov.M.txt",
		log2ratio="results_processed/M_F_norm_coverage_ratio_log2.bed.txt"
	threads: 2
	shell:
		"""
		seqtk comp {input.fa} | awk '{{print $1"\t"$2}}' > genomefile.txt
		bedtools makewindows -w {windowsize} -g genomefile.txt > windows.bed
		rm genomefile.txt

		# separate calls for M and F
		cat {input.popmap} | awk 'NF {{if($2==1) print $1}}' | tr "\\n" "\\t" | awk 'BEGIN {{ OFS="\\t"}}; {{$1=$1; print "chr","start","stop",$0}}' > {output.mcov}
		cat {input.popmap} | awk 'NF {{if($2==2) print $1}}' | tr "\\n" "\\t" | awk 'BEGIN {{ OFS="\\t"}}; {{$1=$1; print "chr","start","stop",$0}}' > {output.fcov}
		( bedtools multicov -bams $( cat {input.popmap} | awk 'BEGIN {{ ORS=" "}}; NF {{if($2==1) print "mapped_reads/"$1".sorted.bam"}}' ) -bed windows.bed >> {output.mcov} )&
		( bedtools multicov -bams $( cat {input.popmap} | awk 'BEGIN {{ ORS=" "}}; NF {{if($2==2) print "mapped_reads/"$1".sorted.bam"}}' ) -bed windows.bed >> {output.fcov} )&
		wait
		rm windows.bed

		################ get M-F norm read coverage ratio log2 (we add 1e-6 to the read counts to avoid zero division)
		awk '{{ for(i=4; i<=NF;i++) j+=$i; print j; j=0 }}' {output.mcov} > total.M
		totm=$( awk '{{sum+=$1;}} END{{print sum;}}' total.M )
		awk -v totalsum="$totm" '{{print ($1*1000000)/totalsum}}' total.M > total.norm.M

		awk '{{ for(i=4; i<=NF;i++) j+=$i; print j; j=0 }}' {output.fcov} > total.F
		totf=$( awk '{{sum+=$1;}} END{{print sum;}}' total.F )
		awk -v totalsum="$totf" '{{print ($1*1000000)/totalsum}}' total.F > total.norm.F

		paste total.norm.M total.norm.F | awk '{{print log((($1+1e-6)/($2+1e-6)))/log(2)}}' > log2_ratio

		paste <(cut -f1,2,3 {output.fcov} ) log2_ratio | tail -n +2 > {output.log2ratio}

		rm total.M total.F total.norm.M total.norm.F log2_ratio
		"""

rule pi_rawstats_M:
	input:
		popmap={popmapfile},
		gzvcf="results_raw/all.post_filter.vcf.gz"
	output:
		"results_raw/pi_raw.M.txt.gz"
	shell:
		"""
		cat {input.popmap} | awk '{{if($2==1) print $1}}' > mpop_pi
		vcftools --gzvcf {input.gzvcf} --keep mpop_pi --max-missing 0.01 --recode --stdout | python scripts/make_denominator_and_numerator_for_pi.py | gzip -c > {output}
		rm mpop_pi
		"""

rule pi_rawstats_F:
	input:
		popmap={popmapfile},
		gzvcf="results_raw/all.post_filter.vcf.gz"
	output:
		"results_raw/pi_raw.F.txt.gz"
	shell:
		"""
		cat {input.popmap} | awk '{{if($2==2) print $1}}' > fpop_pi
		vcftools --gzvcf {input.gzvcf} --keep fpop_pi --max-missing 0.01 --recode --stdout | python scripts/make_denominator_and_numerator_for_pi.py | gzip -c > {output}
		rm fpop_pi
		"""

rule pi_windowed:
	input:
		pim="results_raw/pi_raw.M.txt.gz",
		pif="results_raw/pi_raw.F.txt.gz",
		fa=genomefile
	output:
		pi_f_bed="results_processed/F.pi.bed.txt",
		pi_m_bed="results_processed/M.pi.bed.txt"
	threads: 4
	shell:
		"""
		# unzip raw
		(gunzip -c {input.pim} > tmp.pi_raw.M.txt )&
		(gunzip -c {input.pif} > tmp.pi_raw.F.txt )&
		wait

		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.pi.txt
		bedtools makewindows -w {windowsize} -g genomefile.pi.txt > windows.pi.bed

		# ugly construct to handle the "bug" that chromosomes in the VCF are not sorted in
		# same order as in the fasta genome file; its unclear how this can happen but it DOES
		cut -f1 tmp.pi_raw.M.txt | uniq | awk '{{print $1"\\t0\\t9999999999999999"}}' > pi_sort_order_in_scores.txt
		# find chroms NOT YET in the list
		comm -13 <(cut -f1 pi_sort_order_in_scores.txt | sort | uniq ) <(cut -f1 genomefile.pi.txt | sort | uniq) | awk '{{print $1"\\t0\\t9999999999999999"}}' >> pi_sort_order_in_scores.txt

		bedtools sort -g pi_sort_order_in_scores.txt -i windows.pi.bed > windows.pi.sort_order_in_scores.bed

		# MALE
		(bedtools map -a windows.pi.sort_order_in_scores.bed -b tmp.pi_raw.M.txt -c 4 -o mean > pi_num_M )&
		(bedtools map -a windows.pi.sort_order_in_scores.bed -b tmp.pi_raw.M.txt -c 5 -o mean > pi_denom_M )&

		# FEMALE
		(bedtools map -a windows.pi.sort_order_in_scores.bed -b tmp.pi_raw.F.txt -c 4 -o mean > pi_num_F )&
		(bedtools map -a windows.pi.sort_order_in_scores.bed -b tmp.pi_raw.F.txt -c 5 -o mean > pi_denom_F )&

		wait

		paste pi_num_M pi_denom_M > pi_both_M
		awk '{{ if($8>0) print $1"\\t"$2"\\t"$3"\\t"$4/$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' pi_both_M > pi_tmp
		bedtools sort -g genomefile.pi.txt -i pi_tmp > {output.pi_m_bed}

		paste pi_num_F pi_denom_F > pi_both_F
		awk '{{ if($8>0) print $1"\\t"$2"\\t"$3"\\t"$4/$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' pi_both_F > pi_tmp
		bedtools sort -g genomefile.pi.txt -i pi_tmp > {output.pi_f_bed}

		rm genomefile.pi.txt windows.pi.bed pi_num_M pi_denom_M pi_num_F pi_denom_F pi_both_M pi_both_F pi_tmp windows.pi.sort_order_in_scores.bed pi_sort_order_in_scores.txt tmp.pi_raw.M.txt tmp.pi_raw.F.txt
		"""

rule dxy_rawstats:
	input:
		popmap={popmapfile},
		gzvcf="results_raw/all.post_filter.vcf.gz"
	output:
		"results_raw/dxy_raw.txt.gz"
	shell:
		"""
		zcat {input.gzvcf} | python scripts/make_denominator_and_numerator_for_dxy.py --popmap {input.popmap} | gzip -c > {output}
		"""

rule dxy_windowed:
	input:
		raw="results_raw/dxy_raw.txt.gz",
		fa=genomefile
	output:
		"results_processed/MvF.dxy.bed.txt"
	threads: 2
	shell:
		"""
		# unzip raw file
		gunzip -c {input.raw} > tmp.dxy_raw.txt

		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.dxy.txt
		bedtools makewindows -w {windowsize} -g genomefile.dxy.txt > windows.dxy.bed

		# ugly construct to handle the "bug" that chromosomes in the VCF are not sorted in
		# same order as in the fasta genome file; its unclear how this can happen but it DOES
		cut -f1 tmp.dxy_raw.txt | uniq | awk '{{print $1"\\t0\\t9999999999999999"}}' > dxy_sort_order_in_scores.txt
		# find chroms NOT YET in the list
		comm -13 <(cut -f1 dxy_sort_order_in_scores.txt | sort | uniq ) <(cut -f1 genomefile.dxy.txt | sort | uniq) | awk '{{print $1"\\t0\\t9999999999999999"}}' >> dxy_sort_order_in_scores.txt

		bedtools sort -g dxy_sort_order_in_scores.txt -i windows.dxy.bed > windows.dxy.sort_order_in_scores.bed

		# correct
		( bedtools map -a windows.dxy.sort_order_in_scores.bed -b tmp.dxy_raw.txt -c 4 -o mean > dxy_num )&
		( bedtools map -a windows.dxy.sort_order_in_scores.bed -b tmp.dxy_raw.txt -c 5 -o mean > dxy_denom )&
		wait

		paste dxy_num dxy_denom > dxy_both

		awk '{{ if($8>0) print $1"\\t"$2"\\t"$3"\\t"$4/$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' dxy_both > dxy_tmp

		bedtools sort -g genomefile.dxy.txt -i dxy_tmp > {output}

		rm genomefile.dxy.txt windows.dxy.bed dxy_num dxy_denom dxy_both dxy_tmp windows.dxy.sort_order_in_scores.bed dxy_sort_order_in_scores.txt tmp.dxy_raw.txt
		"""


rule netdiv_windowed:
	input:
		dxy="results_processed/MvF.dxy.bed.txt",
		pi_f_bed="results_processed/F.pi.bed.txt",
		pi_m_bed="results_processed/M.pi.bed.txt",
		fa=genomefile
	output:
		netdiv_M="results_processed/MvF.netdiv_M.bed.txt",
		netdiv_F="results_processed/MvF.netdiv_F.bed.txt"
	shell:
		"""
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.netdiv.txt
		bedtools makewindows -w {windowsize} -g genomefile.netdiv.txt > windows.netdiv.bed

		paste {input.dxy} {input.pi_m_bed} > netdiv_prep
		awk '{{ if($8!="NA" && $4!="NA") print $1"\\t"$2"\\t"$3"\\t"$4-$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' netdiv_prep > {output.netdiv_M}

		paste {input.dxy} {input.pi_f_bed} > netdiv_prep
		awk '{{ if($8!="NA" && $4!="NA") print $1"\\t"$2"\\t"$3"\\t"$4-$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' netdiv_prep > {output.netdiv_F}

		rm genomefile.netdiv.txt windows.netdiv.bed netdiv_prep
		"""


rule fst_rawstats:
	input:
		popmap={popmapfile},
		gzvcf="results_raw/variant_sites.vcf.gz"
	output:
		"results_raw/fst_raw.txt.gz"
	shell:
		"""
		cat {input.popmap} | awk '{{if($2==1) print $1}}' > mpop_fst
		cat {input.popmap} | awk '{{if($2==2) print $1}}' > fpop_fst

		# exclude singletons..https://www.biostars.org/p/278646/
		# also drop any "-nan" outputs immediately
		vcftools --gzvcf {input.gzvcf} --mac 2 --weir-fst-pop mpop_fst --weir-fst-pop fpop_fst --stdout | awk '{{print $1"\\t"$2"\\t"$2"\\t"$3}}' | awk '{{if ( !($4=="-nan") ) print}}' | gzip -c > {output}
		rm mpop_fst fpop_fst
		"""

rule fst_windowed:
	input:
		fraw="results_raw/fst_raw.txt.gz",
		fa=genomefile
	output:
		"results_processed/MvF.fst.bed.txt"
	shell:
		"""
		## we get an "average of ratios" rather than the "ratio of averages"..
		## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3759727/

		# unzip raw
		gunzip -c {input.fraw} > tmp.fst_raw.txt

		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.fst.txt
		bedtools makewindows -w {windowsize} -g genomefile.fst.txt > windows.fst.bed

		fstlines=$(cat tmp.fst_raw.txt | wc -l)
		if [ "$fstlines" -gt 2 ]; then
			bedtools map -a windows.fst.bed -b tmp.fst_raw.txt -g genomefile.fst.txt -c 4 -o mean > {output}
		else
			# fst file was empty, because all values in genome were NA. Make an NA windowed-file to allow downstream to proceed anyway.
			cat windows.fst.bed | awk '{{print $0"\\tNA"}}' > {output}
		fi

		rm genomefile.fst.txt windows.fst.bed tmp.fst_raw.txt
		"""


rule GWAS_plink_raw:
	input:
		popmap={popmapfile},
		gzvcf="results_raw/variant_sites.vcf.gz",
		fa=genomefile
	output:
		"results_raw/plink.assoc_results.assoc.raw.txt.gz"
	resources:
		mem_mb=20000,
		cpus=4
	shell:
		"""
		# PLINK has a few problematic behaviours:
		# - VCF gets converted to temporary plink format files with hard-coded names.. => RACE CONDITON when multiple instances of plink run
		# - it will by default use ALL-1 CPUs => must use argument --threads
		# - By default, PLINK 1.9 tries to reserve half of your system's RAM for its main workspace. => --memory <main workspace size, in MB>

		DIR="tmpdir_plink_GWAS"
		if [ -d "$DIR" ]; then
		  rm -r $DIR
		fi
		mkdir -p tmpdir_plink_GWAS
		cd tmpdir_plink_GWAS

		cat ../{input.popmap} | awk 'NF {{print $1"\\t"$1"\\t"$2}}' > thepheno
		plink --vcf ../{input.gzvcf} --double-id --allow-extra-chr --pheno thepheno --allow-no-sex --assoc --out assoc_results	--threads {resources.cpus} --memory {resources.mem_mb}

		cat assoc_results.assoc | gzip -c > ../{output}
		cd ../
		rm -r tmpdir_plink_GWAS
		"""


rule GWAS_plink_windows:
	input:
		fa=genomefile,
		rawgwas="results_raw/plink.assoc_results.assoc.raw.txt.gz"
	output:
		"results_processed/plink.assoc_results.significant.window.bed.txt",
	shell:
		"""
		zcat {input.rawgwas} | tr -s " " "\\t" > plink.assoc_results.raw.txt

		# Benjamini-Hochberg FDR procedure
		bash scripts/PtoFDR.sh plink.assoc_results.raw.txt 9 > assoc.fdr.txt
		cat assoc.fdr.txt | grep -v "CHISQ" | awk '{{if($9!="NA") print }}' | awk '{{if($10<=0.1) print}}' > fdr10perc.txt
		cat fdr10perc.txt | cut -f1,3,10 | awk '{{print $1"\\t"$2"\\t"$2}}' > plink.assoc_results.significant.bed

		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.assoc.txt
		bedtools makewindows -w {windowsize} -g genomefile.assoc.txt > windows.assoc.bed
		bedtools coverage -counts -a windows.assoc.bed -b plink.assoc_results.significant.bed > {output}

		rm plink.assoc_results.raw.txt assoc.fdr.txt fdr10perc.txt plink.assoc_results.significant.bed genomefile.assoc.txt windows.assoc.bed
		"""

rule LD_plink_raw:
	input:
		gzvcf="results_raw/variant_sites.vcf.gz",
		fa=genomefile
	output:
		"results_raw/ld_clean.sorted.bed.gz"
#		"results_processed/LD_r2_averaged_per_window.bed.txt"
	resources:
		mem_mb=20000,
		cpus=4
	shell:
		"""
		# PLINK has a few problematic behaviours:
		# - VCF gets converted to temporary plink format files with hard-coded names.. => RACE CONDITON when multiple instances of plink run
		# - it will by default use ALL-1 CPUs => must use argument --threads
		# - By default, PLINK 1.9 tries to reserve half of your system's RAM for its main workspace. => --memory <main workspace size, in MB>

		DIR="tmpdir_plink_LD"
		if [ -d "$DIR" ]; then
		  rm -r $DIR
		fi
		mkdir -p tmpdir_plink_LD
		cd tmpdir_plink_LD

		# apply a MAF filter 0.2 before calculating LD as rare alleles imply very high LD by definition; they are thus useless for our purpose.
		# also, thin variants to at least 1000 bp distance; (this in combination with ld-window 500 resulted in best "LD island" visibility in Leucadendron rubrum data)
		vcftools --gzvcf ../{input.gzvcf} --thin 1000 --maf 0.2 --recode --stdout | bgzip -c > forplinkld.vcf.gz
		# vcftools --gzvcf {input.gzvcf} --maf 0.1 --recode --stdout > forplinkld.vcf
		# plink --vcf forplinkld.vcf --double-id --allow-extra-chr --r2 --ld-window-r2 0.0 --ld-window 5 --ld-window-kb 2
		# rm forplinkld.vcf

		plink --vcf forplinkld.vcf.gz --double-id --allow-extra-chr --maf 0.2 --r2 --ld-window-r2 0.0 --ld-window 500 --threads {resources.cpus} --memory {resources.mem_mb}

		# --ld-window X	compute LD only for pairs that are at most X SNPs apart (default 10)
		# --ld-window-kb X	compute LD only for pairs that are at most X kb apart (default 1000 kb)
		# --ld-window-r2 X	minimum r2 to report, else omit from output. Default = 0.2

		#cat plink.ld | tr -s ' ' '\\t' | cut -f1,2,4,5,7 | tail -n +2 | awk '{{ print $1"\\t"$2"\\t"$2"\\t"$5"\\n"$3"\\t"$4"\\t"$4"\\t"$5  }}' | sort --parallel {resources.cpus} -S 2G -T . -k1,1 -k2,2n | gzip -c > ../{output}
		# interpolate the position of the LD value as the midpoint between the two variants from which it is calculated.
		cat plink.ld | tr -s ' ' '\\t' | cut -f1,2,4,5,7 | awk '{{if($1==$3) print $1"\\t"int( ($2+$4)/2 )"\\t"int( ($2+$4)/2 )"\\t"$5 }}' > tmp_plink_before_sort
		sort --parallel {resources.cpus} -S 2G -T . -k1,1 -k2,2n tmp_plink_before_sort | gzip -c > ../results_raw/ld_clean.sorted.bed.gz

		cd ../
		rm -r tmpdir_plink_LD
		"""

rule LD_plink_windows:
	input:
		ldraw="results_raw/ld_clean.sorted.bed.gz",
		fa=genomefile
	output:
		"results_processed/LD_r2_averaged_per_window.bed.txt"
	resources:
		mem_mb=20000,
		cpus=4
	shell:
		"""
		# create windows, sort chroms lexicographically (same as the LD output)
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.ld.txt
		bedtools makewindows -w {windowsize} -g genomefile.ld.txt > tmpldwindows
		sort --parallel {resources.cpus} -S 2G -T . -k1,1 -k2,2n tmpldwindows > windows.ld.bed
		rm tmpldwindows

		gunzip --stdout {input.ldraw} > ld_clean.sorted.bed

		# get the mean LD per window
		bedtools map -a windows.ld.bed -b ld_clean.sorted.bed -c 4 -o mean > tmpoutperwindow.txt

		# sort chroms back to the same order as the genome file:
		bedtools sort -g genomefile.ld.txt -i tmpoutperwindow.txt > {output}

		rm ld_clean.sorted.bed genomefile.ld.txt tmpoutperwindow.txt windows.ld.bed
		"""



rule calc_indiv_het:
	input:
		popm={popmapfile},
		gzvcf="results_raw/variant_sites.vcf.gz"
	output:
		mhet="results_raw/het_males.txt.gz",
		fhet="results_raw/het_females.txt.gz"
	shell:
		"""
		# calculate heterozygosity for GLOBAL variants, i.e. M or F may be fixed for a variant but we calculate heterozygosity anyway.
		# do not report het for sites where no genotypes are present.
		vcftools --gzvcf {input.gzvcf} --keep <( cat {input.popm} | awk '{{if($2==1) print $1}}' ) --hardy --stdout | awk -F'[\\t/]' 'NR>1 {{if(($3+$4+$5)>0) print $1"\\t"$2"\\t"$2"\\t"$4/($3+$4+$5) }}' | gzip -c > {output.mhet}
		vcftools --gzvcf {input.gzvcf} --keep <( cat {input.popm} | awk '{{if($2==2) print $1}}' ) --hardy --stdout | awk -F'[\\t/]' 'NR>1 {{if(($3+$4+$5)>0) print $1"\\t"$2"\\t"$2"\\t"$4/($3+$4+$5) }}' | gzip -c > {output.fhet}
		"""

rule indiv_het_per_windows:
	input:
		fa=genomefile,
		mhet="results_raw/het_males.txt.gz",
		fhet="results_raw/het_females.txt.gz"
	output:
		mhetwindows="results_processed/het_males.bed.txt",
		fhetwindows="results_processed/het_females.bed.txt"
	shell:
		"""
		# unzip raw
		(gunzip -c {input.mhet} > tmp.het_males.txt)&
		(gunzip -c {input.fhet} > tmp.het_females.txt)&
		wait

		# create windows, sort chroms lexicographically (same as the LD output)
		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.het.txt
		bedtools makewindows -w {windowsize} -g genomefile.het.txt > windows.het.bed

		# get the mean HET per window
		bedtools map -a windows.het.bed -b tmp.het_males.txt -c 4 -o mean > {output.mhetwindows}
		bedtools map -a windows.het.bed -b tmp.het_females.txt -c 4 -o mean > {output.fhetwindows}

		rm genomefile.het.txt windows.het.bed tmp.het_males.txt tmp.het_females.txt
		"""

rule on:
	input:
		gzvcf="results_raw/all.post_filter.vcf.gz",
		popm={popmapfile}
	output:
		"results_raw/gametolog_candidate_alleles.allsites.vcf.gz"
	shell:
		"""
		python scripts/pseudo_phase_gametologs.py {input.popm} {input.gzvcf} | bgzip -c > {output}
		"""

rule score_XY_gametolog_divergence:
	input:
		"results_raw/gametolog_candidate_alleles.allsites.vcf.gz"
	output:
		"results_raw/gametolog_candidate_alleles.XY_divergence.rawstats.txt.gz"
	shell:
		"""
		echo -e "XY_Y_like\\tXY_Y_like" > xypopm
		echo -e "XY_X_like\\tXY_X_like" >> xypopm

		echo -e "XY_X_like" > xypop
		echo -e "XY_Y_like" >> xypop

		vcftools --gzvcf {input} --keep xypop --recode --stdout | python scripts/make_denominator_and_numerator_for_dxy.py --popmap xypopm | gzip -c > {output}
		rm xypopm xypop
		"""

rule score_ZW_gametolog_divergence:
	input:
		"results_raw/gametolog_candidate_alleles.allsites.vcf.gz"
	output:
		"results_raw/gametolog_candidate_alleles.ZW_divergence.rawstats.txt.gz"
	shell:
		"""
		echo -e "ZW_W_like\\tZW_W_like" > zwpopm
		echo -e "ZW_Z_like\\tZW_Z_like" >> zwpopm

		echo -e "ZW_W_like" > zwpop
		echo -e "ZW_Z_like" >> zwpop

		vcftools --gzvcf {input} --keep zwpop --recode --stdout | python scripts/make_denominator_and_numerator_for_dxy.py --popmap zwpopm | gzip -c > {output}
		rm zwpopm zwpop
		"""

rule XY_div_windows:
	input:
		raw="results_raw/gametolog_candidate_alleles.XY_divergence.rawstats.txt.gz",
		fa=genomefile
	output:
		"results_processed/gametolog_candidate_alleles.XY_divergence.windows.bed.txt"
	threads: 2
	shell:
		"""
		# unzip raw
		gunzip -c {input.raw} > tmp.gametolog_candidate_alleles.XY_divergence.rawstats.txt

		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.xygametologs.txt
		bedtools makewindows -w {windowsize} -g genomefile.xygametologs.txt > windows.xygametologs.bed

		(bedtools map -a windows.xygametologs.bed -b tmp.gametolog_candidate_alleles.XY_divergence.rawstats.txt -g genomefile.xygametologs.txt -c 4 -o mean > XY_dxy_num )&
		(bedtools map -a windows.xygametologs.bed -b tmp.gametolog_candidate_alleles.XY_divergence.rawstats.txt -g genomefile.xygametologs.txt -c 5 -o mean > XY_dxy_denom )&
		wait

		paste XY_dxy_num XY_dxy_denom > XY_dxy_both

		awk '{{ if($8>0) print $1"\\t"$2"\\t"$3"\\t"$4/$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' XY_dxy_both > {output}

		rm genomefile.xygametologs.txt windows.xygametologs.bed XY_dxy_denom XY_dxy_num XY_dxy_both	tmp.gametolog_candidate_alleles.XY_divergence.rawstats.txt
		"""

rule ZW_div_windows:
	input:
		raw="results_raw/gametolog_candidate_alleles.ZW_divergence.rawstats.txt.gz",
		fa=genomefile
	output:
		"results_processed/gametolog_candidate_alleles.ZW_divergence.windows.bed.txt"
	threads: 2
	shell:
		"""
		# unzip raw
		gunzip -c {input.raw} > tmp.gametolog_candidate_alleles.ZW_divergence.rawstats.txt

		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.zwgametologs.txt
		bedtools makewindows -w {windowsize} -g genomefile.zwgametologs.txt > windows.zwgametologs.bed

		( bedtools map -a windows.zwgametologs.bed -b tmp.gametolog_candidate_alleles.ZW_divergence.rawstats.txt -g genomefile.zwgametologs.txt -c 4 -o mean > ZW_dxy_num )&
		( bedtools map -a windows.zwgametologs.bed -b tmp.gametolog_candidate_alleles.ZW_divergence.rawstats.txt -g genomefile.zwgametologs.txt -c 5 -o mean > ZW_dxy_denom )&
		wait

		paste ZW_dxy_num ZW_dxy_denom > ZW_dxy_both

		awk '{{ if($8>0) print $1"\\t"$2"\\t"$3"\\t"$4/$8 ; else print $1"\\t"$2"\\t"$3"\\tNA" }}' ZW_dxy_both > {output}

		rm genomefile.zwgametologs.txt windows.zwgametologs.bed ZW_dxy_denom ZW_dxy_num ZW_dxy_both	tmp.gametolog_candidate_alleles.ZW_divergence.rawstats.txt
		"""


rule plot_all:
	input:
		regions=regions_for_plot_bed,
		a="results_processed/M_F_norm_coverage_ratio_log2.bed.txt",
		b="results_processed/kmers.F_specific.per_window.bed.txt",
		c="results_processed/kmers.M_specific.per_window.bed.txt",
		d="results_processed/F.pi.bed.txt",
		e="results_processed/M.pi.bed.txt",
		f="results_processed/MvF.fst.bed.txt",
		g="results_processed/MvF.dxy.bed.txt",
		h="results_processed/LD_r2_averaged_per_window.bed.txt",
		i="results_processed/plink.assoc_results.significant.window.bed.txt",
		j="results_processed/gametolog_candidate_alleles.ZW_divergence.windows.bed.txt",
		k="results_processed/gametolog_candidate_alleles.XY_divergence.windows.bed.txt",
		l="results_processed/het_males.bed.txt",
		m="results_processed/het_females.bed.txt",
		n="results_processed/MvF.netdiv_F.bed.txt",
		o="results_processed/MvF.netdiv_M.bed.txt"
	output:
		stats="results_processed/allstats.txt",
		regionstats="results_processed/region.stats.txt",
		regionplot="results_processed/region.stats.plots.pdf"
	shell:
		"""
		cat {input.a} > {output.stats}

		for i in {input.b} {input.c} {input.d} {input.e} {input.f} {input.g} {input.h} {input.i} {input.j} {input.k} {input.l} {input.m} {input.n} {input.o} ; do
			paste {output.stats} <(cut -f4 $i ) > tmpst
			mv tmpst {output.stats}
		done

		bedtools intersect -wa -a {output.stats} -b {input.regions} > {output.regionstats}

		Rscript scripts/plot_allstats.R {output.regionstats} {output.regionplot}
		"""


rule count_reads_input_and_mapping:
	input:
		indexfile="mapped_reads/{sample}.sorted.bam.bai",
		bamfile="mapped_reads/{sample}.sorted.bam",
		reads=get_fastq_sample_ALL
	output:
		temp( "results_raw/{sample}.raw_and_mapped_reads_report.txt" )
	threads: 3
	run:
		cmd = "samtools stats " + input.bamfile + " > " + input.bamfile + ".stats.txt"
		os.system( cmd )
		mean_insert_size = "NA"
		with open(input.bamfile +".stats.txt", "r") as I:
			for line in I:
				if line.startswith("SN	reads mapped:"):
					total_mapped = int(line.strip("\n").split()[-1])
				if line.startswith("SN	insert size average:"):
					mean_insert_size = line.strip("\n").split()[-1]
		os.system( "rm " +  input.bamfile +".stats.txt")
		# reconstruct the sample name...
		samplename = input.bamfile.split("/")[1].split(".sorted.bam")[0]
		# now count all the reads
		sum_of_read_fastq_lines = 0
		for f in list(input.reads):
			os.system( "pigz -p3 -dc {0} | wc -l > tmp.{1}.readcount".format(f,samplename) )
			with open( "tmp.{}.readcount".format(samplename), "r" ) as an_infile:
				sum_of_read_fastq_lines += int(an_infile.readline().strip())
			os.system( "rm tmp.{}.readcount".format(samplename) )
		n_reads = float(sum_of_read_fastq_lines)/4.0
		mapping_rate = float(total_mapped)/n_reads
		with open(str(output), "w") as O:
			O.write(samplename + "\t" + str(n_reads) + "\t" + str(total_mapped) + "\t" + str(round(mapping_rate,3)) + "\t" + mean_insert_size + "\n")

rule collect_mapping_report:
	input:
		expand("results_raw/{sample}.raw_and_mapped_reads_report.txt", sample=SAMPLES)
	output:
		"results_raw/mapping_statistics_report.txt"
	shell:
		"""
		echo -e "counts are single reads, not pairs of reads, i.e. in PE data fwd and rev are counted separately" > {output}
		echo -e "sample\traw_reads\tmapped_reads\tmapping_rate\tinsert_size_average" >> {output}
		cat results_raw/*.raw_and_mapped_reads_report.txt >> {output}
		"""


rule dump_kmer_fastas_and_match_kmers_to_genome:
	input:
		popm={popmapfile},
		fa="bwa_index/reference.fa",
		gidx="bwa_index/reference.fa.bwt",
		fkmers="results_raw/F_specific_kmer.kmc_suf",
		mkmers="results_raw/M_specific_kmer.kmc_suf"
	output:
		fsp="results_processed/kmers.F_specific.per_window.bed.txt",
		msp="results_processed/kmers.M_specific.per_window.bed.txt"
	threads: 12
	shell:
		"""
		# transform kmer databases to fasta
		kmc_tools transform results_raw/F_specific_kmer dump /dev/stdout | awk '{{print ">mer\\n"$1}}' | pigz -p{threads} -c > results_raw/F_specific_kmer.fa.gz
		kmc_tools transform results_raw/M_specific_kmer dump /dev/stdout | awk '{{print ">mer\\n"$1}}' | pigz -p{threads} -c > results_raw/M_specific_kmer.fa.gz

		# remove the kmer databases, keep only the fasta dump:
		rm results_raw/F_specific_kmer.kmc_*
		rm results_raw/M_specific_kmer.kmc_*

		# searching for imperfect matches of 21-mers: https://bioinformatics.stackexchange.com/questions/7298/mapping-heteryzygous-kmers-on-a-genome
		# The -k 21 says to use a minimum seed length of 21, which forces exact matches. The -T 21 requires a minimum score of 21, which also enforces an exact match. The -a parameter reports all matches, since a "best" match doesn't make sense in this situation. The -c parameter limits how many matches are reported, which may need to be adjusted depending on how repetitive the 21-mer is.
		# we match IMPERFECTLY = allowing 1 mismatch in 21 bp kmers. => 4.76% divergence

		# use samtools view -F 4 : removes unmapped

		bwa mem -t {threads} -k 20 -T 20 -a -c 5000 {input.fa} results_raw/F_specific_kmer.fa.gz | samtools view -F 4 -b - > F_specific_kmer.bam
		samtools sort -T F_specific_kmer -O bam F_specific_kmer.bam > F_specific_kmer.sorted.bam
		samtools index F_specific_kmer.sorted.bam
		rm F_specific_kmer.bam #results_raw/F_specific_kmer.fa.gz

		bwa mem -t {threads} -k 20 -T 20 -a -c 5000 {input.fa} results_raw/M_specific_kmer.fa.gz | samtools view -F 4 -b - > M_specific_kmer.bam
		samtools sort -T M_specific_kmer -O bam M_specific_kmer.bam > M_specific_kmer.sorted.bam
		samtools index M_specific_kmer.sorted.bam
		rm M_specific_kmer.bam #results_raw/M_specific_kmer.fa.gz

		seqtk comp {input.fa} | awk '{{print $1"\\t"$2}}' > genomefile.kmer.txt
		bedtools makewindows -w {windowsize} -g genomefile.kmer.txt > windows.kmer.bed
		rm genomefile.kmer.txt
		bedtools multicov -bams F_specific_kmer.sorted.bam -bed windows.kmer.bed > {output.fsp}
		bedtools multicov -bams M_specific_kmer.sorted.bam -bed windows.kmer.bed > {output.msp}
		rm windows.kmer.bed M_specific_kmer.sorted.bam M_specific_kmer.sorted.bam.bai F_specific_kmer.sorted.bam F_specific_kmer.sorted.bam.bai
		"""


rule KMC_get_group_specifics:
	input:
		f="results_raw/group_female.kmc_suf",
		m="results_raw/group_male.kmc_suf"
	output:
		fkmers="results_raw/F_specific_kmer.kmc_suf",
		mkmers="results_raw/M_specific_kmer.kmc_suf"
	run:
		# What is the N of each group?
		N_males = 0
		N_females = 0
		with open(popmapfile, "r") as I:
			for line in I:
				if len(line) > 1:
					if line.strip().split("\t")[1] == "1":
						N_males += 1
					elif line.strip().split("\t")[1] == "2":
						N_females += 1
		print(N_males, N_females)
		os.system("kmc_tools simple results_raw/group_male -ci"+ str( int(0.8*float(N_males))  ) + " results_raw/group_female -ci"+ str( int(0.2*float(N_females))  ) + " kmers_subtract results_raw/M_specific_kmer")
		os.system("kmc_tools simple results_raw/group_female -ci"+ str( int(0.8*float(N_females))  ) + " results_raw/group_male -ci"+ str( int(0.2*float(N_males))  ) + " kmers_subtract results_raw/F_specific_kmer")
		# DO NOT remove inputs
		# os.system(" rm results_raw/group_male.kmc*  results_raw/group_female.kmc* ")

rule count_kmer_group_level:
	input:
		kmc_dumps=expand("results_raw/{sample}.kmc_dump.fa.gz", sample=SAMPLES),
		popmap={popmapfile}
	output:
		f="results_raw/group_female.kmc_suf",
		m="results_raw/group_male.kmc_suf"
	threads: 24
	shell:
		"""
		# count kmers that occur at least 1 times (-ci), and count them up to 1000 million (-cs)
		cat {input.popmap} | awk '{{if($2==1) print "results_raw/"$1".kmc_dump.fa.gz"}}' >input_file_names.male_fastas
		mkdir -p ./kmc_tempdir.male_fasta_counting
		kmc -m10 -sm -k21 -fa -ci1 -cs1000 -t{threads} @input_file_names.male_fastas results_raw/group_male ./kmc_tempdir.male_fasta_counting
		rm -r ./kmc_tempdir.male_fasta_counting
		rm input_file_names.male_fastas

		cat {input.popmap} | awk '{{if($2==2) print "results_raw/"$1".kmc_dump.fa.gz"}}' >input_file_names.female_fastas
		mkdir -p ./kmc_tempdir.female_fasta_counting
		kmc -m10 -sm -k21 -fa -ci1 -cs1000 -t{threads} @input_file_names.female_fastas results_raw/group_female ./kmc_tempdir.female_fasta_counting
		rm -r ./kmc_tempdir.female_fasta_counting
		rm input_file_names.female_fastas
		"""


rule KMC_count_and_dump_to_fasta:
	input:
		get_fastq_sample_ALL
	output:
		temp( "results_raw/{sample}.kmc_dump.fa.gz" )
	threads: 12
	shell:
		"""
		# count kmers that occur at least 2 times (-ci), and count them up to 2
		echo {input} | tr " " "\n" >input_file_names.{wildcards.sample}
		mkdir -p ./kmc_tempdir.{wildcards.sample}
		kmc -m10 -sm -k21 -fq -ci2 -cs2 -t{threads} @input_file_names.{wildcards.sample} results_raw/kmer_{wildcards.sample} ./kmc_tempdir.{wildcards.sample}
		kmc_tools transform results_raw/kmer_{wildcards.sample} dump /dev/stdout | awk '{{print ">mer\\n"$1}}' | pigz -p{threads} -c > {output}
		rm input_file_names.{wildcards.sample}
		rm -r ./kmc_tempdir.{wildcards.sample}
		rm results_raw/kmer_{wildcards.sample}.kmc_suf results_raw/kmer_{wildcards.sample}.kmc_pre
		"""
