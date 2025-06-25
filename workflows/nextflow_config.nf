// Simple Nextflow configuration for spatial transcriptomics

manifest {
    name = 'spatial-transcriptomics'
    description = 'Spatial transcriptomics analysis pipeline'
    version = '1.0.0'
}

// Default parameters
params {
    input_dir = '/workspace/data/example_data'
    outdir = '/workspace/results'
    count_file = 'Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5'
}

// Process configuration
process {
    // Default resources
    cpus = 2
    memory = '4 GB'
    time = '1h'

    // Error strategy
    errorStrategy = 'retry'
    maxRetries = 2
}

// Execution profiles
profiles {
    standard {
        executor.name = 'local'
        executor.cpus = 4
        executor.memory = '8 GB'
    }

    docker {
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
    }
}

// Reporting
timeline {
    enabled = true
    file = "${params.outdir}/pipeline_info/timeline.html"
}

report {
    enabled = true
    file = "${params.outdir}/pipeline_info/report.html"
}