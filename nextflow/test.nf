#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workDir = "/fs/mwa_base/test"

params.cand = "Blind"
params.fits_file = "./mock/"
params.zmax = 0
params.out_dir = "./mock/out"
params.publish_all_prepfold = true

fits_file = Channel.fromPath( "${params.fits_file}", checkIfExists: true )
fits_file.view( it -> "Running search on ${it}" )


process get_centre_freq {
    output:
    file 'centre_freq.txt'

    """
    echo 'mock by using an existing get_centre_freq.txt'
    """
}

process ddplan {
    input:
    tuple val(name), val(centre_freq)

    output:
    file 'DDplan.txt'

    """
    echo 'ddplan -> mock using an existing DDplan.txt'
    """
}

process search_dd_fft_acc {
    input:
    tuple val(name), val(dm_values), file(fits_files)

    output:
    tuple val(name), file("*ACCEL_${params.zmax}"), file("*.inf"), file("*.subSpS"), file('*.cand')

    """
    echo 'search_dd_fft_acc -> TBD'
    """
}

process accelsift {
    input:
    tuple val(name), file(accel_inf_single_pulse)

    output:
    tuple val(name), file("cands_*greped.txt")

    """
    echo "accelsift -> TBD"
    """
}

process single_pulse_searcher {
    publishDir params.out_dir, mode: 'copy'

    input:
    tuple val(name), file(sps), file(fits)

    output:
    file "*pdf" optional true
    file "*.SpS" optional true

    if ( "$HOSTNAME".startsWith("farnarkle") || "$HOSTNAME".startsWith("x86") ||\
         "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/sps/sps.sif"
    }
    else {
        container = "nickswainston/sps"
    }
    """
     cat *.subSpS > ${name}.SpS
    #-SNR_min 4 -SNR_peak_min 4.5 -DM_cand 1.5 -N_min 3
     single_pulse_searcher.py -fits ${fits} -no_store -plot_name ${name}_sps.pdf ${name}.SpS
    """
}

process prepfold {
    publishDir params.out_dir, mode: 'copy', enabled: params.publish_all_prepfold

    input:
    tuple val(cand_line), file(cand_file), file(cand_inf), file(fits_files)

    output:
    file "*pfd*"

    """
    echo "${cand_line.split()}"
    """
}


workflow pulsar_search {
    take:
        name_fits_files
    main:
        get_centre_freq()
        ddplan( name_fits_files.map{ it -> it[0] }.combine(get_centre_freq.out.splitCsv()) )
        search_dd_fft_acc( // combine the fits files and ddplan with the matching name key (candidateName_obsid_pointing)
                           ddplan.out.splitCsv().map{ it -> [ it[0], [ it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] ] }.\
                           concat(name_fits_files).groupTuple().\
                           // Find for each ddplan match that with the fits files and the name key then change the format to [val(name), val(dm_values), file(fits_files)]
                           map{ it -> [it[1].init(), [[it[0], it[1].last()]]].combinations() }.flatMap().\
                           map{ it -> [it[1][0], it[0], it[1][1]]} )
        // Get all the inf, ACCEL and single pulse files and sort them into groups with the same name key
        accelsift( search_dd_fft_acc.out.map{ it -> [it[0], [it[1]].flatten().findAll { it != null } + \
                                                            [it[2]].flatten().findAll { it != null }] }.\
                   groupTuple( size: total_dm_jobs, remainder: true ).map{ it -> [it[0], it[1].flatten()]} )
        single_pulse_searcher( search_dd_fft_acc.out.map{ it -> [it[0], [it[3]].flatten().findAll { it != null }] }.\
                               groupTuple( size: total_dm_jobs, remainder: true ).map{ it -> [it[0], it[1].flatten()]}.\
                               // Add fits files
                               concat(name_fits_files).groupTuple( size: 2 ).map{ it -> [it[0], it[1][0], it[1][1]]} )
        prepfold( name_fits_files.cross(
                  // Group all the accelsift lines together
                  accelsift.out.map{ it -> it[1] }.splitCsv().flatten().map{ it -> [it.split()[0].split("_ACCEL")[0], it ] }.cross(
                  // Group all the .cand and .inf files by their base names
                  search_dd_fft_acc.out.map{ it -> [it[2]].flatten().findAll { it != null } }.
                  flatten().map{ it -> [it.baseName.split(".inf")[0], it ] }.concat(
                  search_dd_fft_acc.out.map{ it -> [it[4]].flatten().findAll { it != null } }.
                  flatten().map{ it -> [it.baseName.split("_ACCEL")[0], it ] }).groupTuple( size: 2 )
                  // match the cand and inf file with each accelsift line and reoraganise
                  ).map{ it -> [it[0][0].split("_DM")[0], [it[0][1], it[1][1][0], it[1][1][1]]] }
                  // Match with fits files and eogranise to val(cand_line), file(cand_file), file(cand_inf), file(fits_files)
                  ).map{ it -> [it[1][1][0], it[1][1][2], it[1][1][1], it[0][1]] } )
    emit:
        accelsift.out
        prepfold.out
}

workflow {
    pulsar_search( fits_file.toSortedList().map{ it -> [ params.cand + '_' + it[0].getBaseName().split("/")[-1].split("_ch")[0], it ] } )
}
