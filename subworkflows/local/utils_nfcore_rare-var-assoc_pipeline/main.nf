//
// Subworkflow with functionality specific to the psuszyns/rare-var-assoc-nf pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { RSCRIPT_BUILD_PHENOTYPES  } from '../../../modules/local/rscript/build_phenotypes'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process JOIN_CASES_AND_CONTROLS {
    label 'process_low'

    input:
    tuple val(meta), path(cases), path(controls)
    val(options)

    output:
    tuple val(meta), path("all.samples"), emit: output_file

    script:
    def fromMatcher = (options =~ /--replace-char-from\s*['"](.?)['"]/)
    def toMatcher = (options =~ /--replace-char-to\s*['"](.?)['"]/)

    def sedArg = "s///g"
    if (fromMatcher.find() && toMatcher.find()) {
        def fromChar = fromMatcher[0][1]
        def toChar = toMatcher[0][1]
        sedArg = "s/${fromChar}/${toChar}/g"
    }
    """
    echo "sedArg = ${sedArg}"
    cut -f1 ${cases} > all.samples
    cut -f1 ${controls} >> all.samples
    sed -i ${sedArg} all.samples
    """
}

process PHENOTYPE_SAMPLES {
    label 'process_low'

    input:
    tuple val(meta), path(phenotype)

    output:
    tuple val(meta), path("all.samples"), emit: output_file

    script:
    """
    cut -f1 ${phenotype} | tail -n +2 > all.samples
    """
}

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input_vcf         //  string: Path to input samplesheet
    input_controls
    input_cases
    input_phenotype

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //
    if (input_phenotype != null) {
        // input_phenotype may be a comma-separated list of files
        input_phenotype = input_phenotype.tokenize(',')
        Channel.fromPath(input_phenotype, checkIfExists: true)
            .map { file ->
                log.info("Processing phenotype file: ${file.name}")
                def matcher = (file.name =~ /.*_dataset_idx_(\d+)_.*\.phenotype\.txt/)
                if (matcher.matches()) {
                    def dataset_idx = matcher[0][1] // Extract the dataset identifier
                    def meta_id = "${params.project_name}_dataset_idx_${dataset_idx}"
                    tuple([id: meta_id], file)
                } else {
                    log.warn("Phenotype file name does not match expected pattern. Using project_name as id.")
                    tuple( [id: params.project_name], file )
                }
            }
            .set { ch_phenotype }

        PHENOTYPE_SAMPLES (
            ch_phenotype
        )
        ch_all_samples = PHENOTYPE_SAMPLES.out.output_file
    } else {
        if (input_controls != null && file(input_controls).isFile() && input_cases != null && file(input_cases).isFile()) {
            Channel.fromPath(input_controls, checkIfExists: true)
                .map { files ->
                    tuple( [id: params.project_name], files ) // Add meta component
                }
                .set { ch_input_controls }

            Channel.fromPath(input_cases, checkIfExists: true)
                .map { files ->
                    tuple( [id: params.project_name], files ) // Add meta component
                }
                .set { ch_input_cases }
            
            r_script_build_phenotypes_ch = Channel.fromPath(params.rscript_build_phenotypes_path, checkIfExists: true)
            RSCRIPT_BUILD_PHENOTYPES (
                ch_input_controls.join(ch_input_cases, by: 0)
                    .combine(r_script_build_phenotypes_ch),
                Channel.value(params.rscript_build_phenotypes_options)
            )
            ch_phenotype  = RSCRIPT_BUILD_PHENOTYPES.out.phenotype
            ch_versions = ch_versions.mix(RSCRIPT_BUILD_PHENOTYPES.out.versions.first())
            
            JOIN_CASES_AND_CONTROLS (
                ch_input_cases.join(ch_input_controls, by: 0),
                Channel.value(params.rscript_build_phenotypes_options)
            )
            ch_all_samples = JOIN_CASES_AND_CONTROLS.out.output_file

        } else {
            throw new RuntimeException("Please provide either a phenotype file or both controls and cases files")
        }
    }

    // Process VCF files and create Cartesian product with phenotype files
    Channel.fromPath(input_vcf, checkIfExists: true)
        .map { vcf -> tuple(vcf) } // Wrap VCF file in a tuple for combining
        .combine(ch_phenotype) // Cartesian product: each VCF with each phenotype
        .map { vcf, pheno_meta, pheno_file ->
            // Create a new tuple with the same meta.id as the phenotype and the VCF file
            tuple([id: pheno_meta.id], vcf)
        }
        .set { ch_input_vcf }

    emit:
    vcf = ch_input_vcf
    phenotype = ch_phenotype
    all_samples = ch_all_samples
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

