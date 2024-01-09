params.binSize = 1000
params.reference
params.metadata
params.pseudoautosome = 6645000

process make_reference_dictionary {
    input:
    file(reference)

    output:
    file("${reference.baseName}.dict")

    script:
    """
    gatk CreateSequenceDictionary \
        -R ${reference} \
        -O ${reference.baseName}.dict
    """
}

process make_reference_index {
    input:
    file(reference)

    output:
    file("${reference}.fai")

    script:
    """
    samtools faidx ${reference}
    """
}

process make_bins_from_reference {
    input:
    file(reference)
    file(index)
    file(dict)

    output:
    file("intervals.interval_list")
    
    script:
    """
    gatk PreprocessIntervals \
        -R ${reference} \
        --bin-length ${params.binSize} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O intervals.interval_list
    """
}

process count_reads_in_bins {
    cpus 1
    executor 'lsf'
    queue 'normal'
    memory 4.GB
    time 12.h

    publishDir "${params.outputDir}/counts", mode: 'copy'
    
    input:
    tuple val(id), file(bam), file(bins), file(reference), file(referenceFai), file(referenceDict)

    output:
    tuple val(id), file("${id}.counts.tsv")

    script: 
    """
    gatk CollectReadCounts \
        -L ${bins} \
        -R ${reference} \
        -imr OVERLAPPING_ONLY \
        -I ${bam[0]} \
        --format TSV \
        -O ${id}.counts.tsv
    gatk \
    """
}

process clean_readcount_file {
    input:
    tuple val(id), path(counts)

    output:
    path "${id}.cleaned_counts.tsv.gz"

    publishDir "${params.outputDir}/cleaned_counts", mode: 'copy'

    script:
    """
    name=\$(grep -E '^@.*SM:' "${counts}" | sed 's/^@.*SM:\s*//')
    grep -v "^@" "${counts}" |\
      awk -v name="\$name" 'BEGIN { OFS="	" } { if (NR == 1) { print "CHROM", \$2, \$3, \$4, "samplename" } else { print \$0, name } }' |\
      gzip > ${id}.cleaned_counts.tsv.gz
    """
}

process merge_readcounts {
    input:
    path tsvs

    output:
    path "readcount_matrices.tsv.gz"

    publishDir "${params.outputDir}/readcount_matrices", mode: 'copy'

    script:
    """
    merge_readcount_files.R \
      --input_path . \
      --output_file readcount_matrices.tsv.gz
    """
}

process make_gc_bias {
    input:
    path bins
    path reference
    path referenceFai

    output:
    path "gc_content.bed"

    publishDir "${params.outputDir}/gc_content", mode: 'copy'

    script:
    """
    bedtools makewindows -g ${referenceFai} -w ${params.binSize} > ivls.bed
    grep -v "^@" ivls.bed |\
      cut -f1-3 |\
      bedtools nuc -fi ${reference} -bed - > gc_content.bed
    """
}

process compute_logrs {
    cpus 2
    executor 'lsf'
    queue 'normal'
    memory 32.GB
    time 12.h

    input:
    path readcount_matrix
    path metadata

    output:
    path "parquet/dataset", emit: dataset
    path "parquet/*.RDS", emit: correctors

    publishDir "${params.outputDir}/compute_logrs", mode: 'link'

    script:
    """
    1.1_compute_tumour_logrs.R \
        --readcounts ${readcount_matrix} \
        --metadata ${metadata} \
        --outputpath parquet
    """
}

process host_cn_filter {
    executor 'lsf'
    queue 'normal'
    memory 16.GB
    time 4.h

    input:
    path inputs
    path metadata
    path modelfile
    val modeltype

    output:
    path "host_cn_filter/*.RDS", emit: rds
    path "host_cn_filter/*.csv", emit: excluded

    publishDir "${params.outputDir}/host_cn_filter", mode: 'link'

    script:
    """
    1.2_host_cn_filter.R \
        --inputpath "." \
        --metadata ${metadata} \
        --modelfile ${modelfile} \
        --modeltype ${modeltype} \
        --outputpath host_cn_filter
    """
}

process host_cn_filter_neuralnet {
    executor 'lsf'
    queue 'normal'
    memory 16.GB
    time 4.h

    input:
    path inputs
    path metadata
    path modelfile
    val modeltype

    output:
    path "host_cn_filter_neuralnet/*.RDS", emit: rds
    path "host_cn_filter_neuralnet/*.csv", emit: excluded

    publishDir "${params.outputDir}/host_cn_filter", mode: 'link'

    script:
    """
    1.2_host_cn_filter.R \
        --inputpath "." \
        --metadata ${metadata} \
        --modelfile ${modelfile} \
        --modeltype ${modeltype} \
        --outputpath host_cn_filter_neuralnet
    """
}

process segment_samples {
    executor 'lsf'
    queue 'normal'
    memory 16.GB
    time 4.h

    input:
    path cn_filter
    path dataset
    path metadata
    path excluded
    each samplename

    output:
    path "segmented/*.metadata.csv", emit: ploidy_metadata
    path "segmented/*.ploidy_assessment.csv", emit: ploidy_assessment
    path "segmented/*.segmentation.csv", emit: segmentation
    path "segmented/*.filtered_segmentation.csv", emit: filtered_segmentation

    publishDir "${params.outputDir}/segmented", mode: 'link'

    script:
    """
    1.3_segment_sample.R \
        dataset/ \
        ${samplename} \
        ${metadata} \
        segmented \
        --excluded_bins ${cn_filter} \
        --pseudoautosome_end ${params.pseudoautosome} \
        --do_ploidy_assessment
    """
}

process combine_ploidy_estimates {
    input:
    path metadata_files

    output:
    path "metadata_with_ploidy.csv"

    publishDir "${params.outputDir}/ploidy"

    script:
    """
    combine_ploidy_estimates.R \
        --input ${metadata_files} \
        --output metadata_with_ploidy.csv
    """
}

process merge_breakpoints {
    executor 'lsf'
    queue 'normal'
    memory 16.GB
    time 4.h

    input:
    path segmented
    path metadata
    val clone

    output:
    path "merged_breakpoints/*"

    publishDir "${params.outputDir}/merged_breakpoints", mode: 'link'

    script:
    """
    1.4_merge_breakpoints.R \
        . \
        ${metadata} \
        ${clone} \
        merged_breakpoints
    """
}

process filter_breakpoints {
    executor 'lsf'
    queue 'normal'
    memory 16.GB
    time 4.h

    input:
    path dataset
    val samplename
    path segmentation
    path metadata

    output:
    path "filter_breakpoints/*"

    publishDir "${params.outputDir}/filter_breakpoints", mode: 'link'

    script:
    """
    1.5_filter_breakpoints.R \
        dataset \
        ${samplename} \
        ${segmentation} \
        ${metadata} \
        filter_breakpoints
    """

}

workflow {
    // Prepare inputs
    reference = Channel.fromPath("${params.reference}", checkIfExists: true)
    metadata = Channel.fromPath("${params.metadata}", checkIfExists: true)
    samplenames = metadata |
        splitCsv(sep: "\t", header: true) |
        map { it -> it.tumour }
    clones = metadata |
        splitCsv(sep: "\t", header: true) |
        map { it -> it.clone }


    dictPath = params.reference.substring(0, params.reference.lastIndexOf('.')) + ".dict"
    faiPath = params.reference + ".fai"

    if (file(dictPath).exists()) {
        referenceDict = Channel.fromPath(dictPath, checkIfExists: true)
    } else {
        referenceDict = make_reference_dictionary(reference)
    }
    
    if (file(faiPath).exists()) {
        referenceFai = Channel.fromPath(faiPath, checkIfExists: true)
    } else {
        referenceFai = make_reference_index(reference)
    }

    bams = Channel.fromFilePairs("${params.inputDir}/*.bam{,.bai}")
    crams = Channel.fromFilePairs("${params.inputDir}/*.cram{,.crai}")
    inputs = bams.concat(crams)

    // Generate the bins and then count reads;
    bins = make_bins_from_reference(reference, referenceFai, referenceDict)
    gc = make_gc_bias(bins, reference, referenceFai)

    countsInput = inputs
        .combine(bins)
        .combine(reference)
        .combine(referenceFai)
        .combine(referenceDict)
    counts = count_reads_in_bins(countsInput)

    cleaned_counts = clean_readcount_file(counts)

    // Merge the readcounts
    merged = merge_readcounts(cleaned_counts.collect())

    // Compute logR
    logrs = compute_logrs(merged, metadata)

    // Filter host CN
    // host_cn_filtered = host_cn_filter(logrs.dataset, metadata, "${params.modelFile}", "${params.modelType}")
    // host_cn_filtered_nn = host_cn_filter_neuralnet(logrs.dataset, metadata, "${params.neuralnet}", "NN")

    // Segment samples
    // segmented = segment_samples(host_cn_filtered.excluded, logrs.dataset, metadata, reference, samplenames)
    // ploidy_metadata = combine_ploidy_estimates(segmented.ploidy_metadata.collect())

    // merge breakpoints
    // merged = merge_breakpoints(segmented.filtered_segmentation.collect(), ploidy_metadata, params.clone)

    // filter breakpoints
    // filtered = filter_breakpoints(logrs.dataset, samplenames, merged, ploidy_metadata)

}
