def classify_sample(sample, files) {

    def meta = [:]
    meta.is_paired = (files instanceof Collection && files.size() == 2)
    meta.id = sample
	meta.pair_id = (files[0].name ==~ /.*R[12]\.(fastq|fq)(\.gz)?/) ? "R" : ""

    return [meta, files]

    if (meta.is_paired) {
        return [meta, files]
    }

    return [meta, [files]]

}
