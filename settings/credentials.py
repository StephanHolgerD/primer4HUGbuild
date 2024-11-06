params = {
    "version": "hg38",
    "data": {
        "reference": "/mnt/data/HG38/GRCh38_latest_genomic.fna",
        "annotation": "/mnt/data/HG38/hg38-p14_annotation.db",
        "coordinates": "/mnt/data/cdot-0.2.1.refseq.grch37_grch38.json.gz",
        "sequences": "/mnt/data/seqrepo/2021-01-29",
        "variation": {
            "dbSNP": "/mnt/data/HG38/GRCh38_latest_dbSNP_all.vcf.gz"
        },
        "chrom_names": "/mnt/data/HG38/chrom_names_hg38.csv"
    },
    "snv_filter": {
        "min_databases": 2,
        "max_snv_len": 10
    },
    "primers": {
        "check_max_num_candidates": 100,
        "mn_3prime_matches": 15,
        "mx_amplicon_len": 4000,
        "mx_amplicon_n": 1
    },
    "blast": {
        "word_size": 13,
        "mx_evalue": 5,
        "n_cpus": 8,
        "mx_blast_hits": 10000,
        "index": "/mnt/data/redux.fna"
    },
    "n_return": 10,
    "burnin_sanger": 30,
    "binding_site": 50,
    "size_min": 18,
    "size_opt": 20,
    "size_max": 27,
    "size_range_PCR": [350, 600],
    "size_range_qPCR": [80, 150],
    "size_range_mRNA": [80, 600],
    "tm_min": 57.0,
    "tm_opt": 60.0,
    "tm_max": 63.0,
    "GC_min": 20.0,
    "GC_max": 80.0,
    "Ns_max": 0,
    "homopolymer_max_len": 4,
    "3prime_max_GC": 4,
    "3prime_stability": 100,
    "Ns_max": 0,
    "salt_monovalent": 50.0,
    "salt_divalent": 1.5,
    "conc_dNTP": 0.6,
    "conc_DNA": 50.0
}
