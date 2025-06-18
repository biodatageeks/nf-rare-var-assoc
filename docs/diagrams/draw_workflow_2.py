from diagrams import Cluster, Diagram
from diagrams.c4 import Person, Container, Database, System, SystemBoundary, Relationship
from diagrams.aws.compute import ElasticContainerServiceService
from diagrams.aws.analytics import AmazonOpensearchService

graph_attr = {
    "splines": "spline",
}

# direction="LR",
with Diagram("Nextflow gene-level associacion pipeline", direction="LR", graph_attr=graph_attr):

    phenotype_file = ElasticContainerServiceService("Phenotype file")
    
    with Cluster("Reports and outputs"):
        eda_diagrams = AmazonOpensearchService("Input Data Characteristics Plots")
        pc_plot = AmazonOpensearchService("Input Data Characteristics Plots")
        html_report = AmazonOpensearchService("HTML Report with Mahattan and QQ plots")
        csv_reports = AmazonOpensearchService("CSV files with results")
        
    
    with Cluster("pVCF Data Ingestion", direction="TB"):
        fasta_file = ElasticContainerServiceService("FASTA reference genome")
        bam_files = ElasticContainerServiceService("BAM/CRAM files")
        gvcf_files = ElasticContainerServiceService("gVCF files")
        input_pvcf_file = ElasticContainerServiceService("Input pVCF file")
        
        deep_variant = Container(
            name="DeepVariant",
            technology="",
            description="",
        )

        merge_gvcfs = Container(
            name="Merge gvcfs",
            technology="",
            description="",
        )
        
        gl_nexus = Container(
            name="GLNexus",
            technology="",
            description="",
        )

        bam_files >> Relationship("") >> deep_variant
        fasta_file >> Relationship("") >> deep_variant
        
        gvcf_files >> Relationship("") >> merge_gvcfs
        deep_variant >> Relationship("") >> merge_gvcfs
        
        merge_gvcfs >> Relationship("") >> gl_nexus
        
        with Cluster("Alternative Input"):
            final_pvcf_file = ElasticContainerServiceService("Final pVCF file")
            
            gl_nexus >> Relationship("") >> final_pvcf_file
            input_pvcf_file >> Relationship("") >> final_pvcf_file
        
    with Cluster("Data Preprocessing", direction="TB"):
        
        exploratory_data_analysis = Container(
            name="Python script",
            technology="",
            description="Exploratory Data Analysis",
        )
        
        bcftools_quality_filtering = Container(
            name="BCFTools",
            technology="",
            description="Genotype quality, depth of coverage filtering, variant normalization, duplicates removal",
        )
        
        vep = Container(
            name="VEP",
            technology="",
            description="Variant annotation",
        )
        
        extract_gnomad_allele_freq = Container(
            name="BCFTools",
            technology="",
            description="Extract gnomad allele frequency data",
        )
        
        plink_missingness_per_pheno = Container(
            name="Plink",
            technology="",
            description="Impute sex, missingness filtering per phenotype",
        )
        
        final_pvcf_file >> Relationship("") >> exploratory_data_analysis
        exploratory_data_analysis >> Relationship("") >> eda_diagrams
        
        final_pvcf_file >> Relationship("") >> bcftools_quality_filtering
        
        bcftools_quality_filtering >> Relationship("") >> vep
        
        vep >> Relationship("") >> extract_gnomad_allele_freq
        vep >> Relationship("") >> plink_missingness_per_pheno

    with Cluster("Data Filtering for Regenie Step 1"):
        
        plink_3 = Container(
            name="Plink",
            technology="",
            description="Minor allele frequency and count filtering, Hardy-Weinberg equilibrium deviations filtering",
        )
        
        with Cluster("Inbreeding Coefficient Filtering", direction="TB"):
            plink_indep_pairwise = Container(
                name="Plink",
                technology="",
                description="--indep-pairwise linkage equilibrium calculation",
            )
            
            plink_het = Container(
                name="Plink",
                technology="",
                description="--het F-coefficient calculation",
            )
            
            plink_f_filtering = Container(
                name="Plink",
                technology="",
                description="Calculate inbreeding outliers in a Python script and then remove them using Plink",
            )
            
            plink_indep_pairwise >> Relationship("") >> plink_het
            plink_het >> Relationship("") >> plink_f_filtering
        
        with Cluster("Principal Component Analysis", direction="TB"):
            
            plink_remove_hild = Container(
                name="Plink",
                technology="",
                description="Remove high-LD regions",
            )
            
            plink_kinship_filtering = Container(
                name="Plink",
                technology="",
                description="Remove related individuals (--king-cutoff)",
            )
            
            plink_pca = Container(
                name="Plink",
                technology="",
                description="Perform PCA",
            )
            
            draw_pc_lot = Container(
                name="Python script",
                technology="",
                description="Draw Principal Components plot",
            )
            
            plink_remove_hild >> Relationship("") >> plink_kinship_filtering
            
            extract_gnomad_allele_freq >> Relationship("") >> plink_pca
            plink_kinship_filtering >> Relationship("") >> plink_pca
            
            plink_pca >> Relationship("") >> draw_pc_lot
            
            draw_pc_lot >> Relationship("") >> pc_plot
        
        
        plink_missingness_per_pheno >> Relationship("") >> plink_3
        plink_3 >> Relationship("") >> plink_indep_pairwise
        plink_f_filtering >> Relationship("") >> plink_remove_hild
        
    with Cluster("Regenie"):
        
        regenie_step_1 = Container(
            name="Regenie",
            technology="",
            description="Step 1",
        )
        
        
        rscript_prepare_annotations = Container(
            name="RScript",
            technology="",
            description="Assign masks",
        )
        
        plink_pgen = Container(
            name="Plink",
            technology="",
            description="Make PGEN/PVAR/PSAM files for optimal Regenie performance",
        )
    
        regenie_step_2 = Container(
            name="Regenie",
            technology="",
            description="Step 2",
        )
    
        phenotype_file >> Relationship("") >> regenie_step_1
        plink_f_filtering >> Relationship("") >> regenie_step_1
        plink_pca >> Relationship("covariates file") >> regenie_step_1
        
        vep >> Relationship("") >> rscript_prepare_annotations
        plink_missingness_per_pheno >> Relationship("") >> plink_pgen
        
        rscript_prepare_annotations >> Relationship("") >> regenie_step_2
        plink_pgen >> Relationship("") >> regenie_step_2
        
    build_reports = Container(
        name="RScript",
        technology="",
        description="Build reports",
    )
    
    regenie_step_2 >> Relationship("") >> build_reports
    eda_diagrams >> Relationship("") >> build_reports
    pc_plot >> Relationship("") >> build_reports
    
    build_reports >> Relationship("") >> html_report
    build_reports >> Relationship("") >> csv_reports
    
    
