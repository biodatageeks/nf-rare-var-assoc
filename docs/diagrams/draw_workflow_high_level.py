from diagrams import Cluster, Diagram
from diagrams.c4 import Person, Container, Database, System, SystemBoundary, Relationship
from diagrams.aws.compute import ElasticContainerServiceService
from diagrams.aws.analytics import AmazonOpensearchService

graph_attr = {
    "splines": "spline",
}

# direction="LR",
with Diagram("Nextflow gene-level associacion pipeline - high level", direction="LR", graph_attr=graph_attr):

    with Cluster("Inputs"):
        phenotype_file = ElasticContainerServiceService("Phenotype file")
        fasta_file = ElasticContainerServiceService("FASTA reference genome")
        bam_files = ElasticContainerServiceService("BAM/CRAM files")
        gvcf_files = ElasticContainerServiceService("gVCF files")
        input_pvcf_file = ElasticContainerServiceService("Input pVCF file")
    
    with Cluster("Reports and outputs"):
        html_report = AmazonOpensearchService("HTML Report with Mahattan \nand QQ plots")
        csv_reports = AmazonOpensearchService("CSV files with results")
        
    data_ingestion = Container(
        name="pVCF Data Ingestion",
        technology="",
        description="",
    )
    
    data_preprocessing = Container(
        name="Data Preprocessing",
        technology="",
        description="",
    )
    
    with Cluster("Data Filtering for Regenie Step 1"):
        
        inbreeding_coefficient_filtering = Container(
            name="Inbreeding Coefficient \nFiltering",
            technology="",
            description="                                            ",
        )
    
        pca = Container(
            name="Principal Component \nAnalysis",
            technology="",
            description="                                            ",
        )
    
    with Cluster("Regenie"):
        
        regenie_step_1 = Container(
            name="Regenie",
            technology="",
            description="Step 1",
        )
    
        regenie_step_2 = Container(
            name="Regenie",
            technology="",
            description="Step 2",
        )
    
    build_reports = Container(
        name="RScript",
        technology="",
        description="Build reports",
    )

    fasta_file >> Relationship("") >> data_ingestion
    bam_files >> Relationship("") >> data_ingestion
    gvcf_files >> Relationship("") >> data_ingestion
    input_pvcf_file >> Relationship("") >> data_ingestion
    phenotype_file >> Relationship("") >> data_preprocessing
    data_ingestion >> Relationship("") >> data_preprocessing
    data_preprocessing >> Relationship("") >> inbreeding_coefficient_filtering
    inbreeding_coefficient_filtering >> Relationship("") >> pca
    inbreeding_coefficient_filtering >> Relationship("") >> regenie_step_1
    data_preprocessing >> Relationship("") >> regenie_step_2
    pca >> Relationship("covariates file") >> regenie_step_1
    pca >> Relationship("covariates file") >> regenie_step_2
    regenie_step_2 >> Relationship("") >> build_reports
    build_reports >> Relationship("") >> html_report
    build_reports >> Relationship("") >> csv_reports

