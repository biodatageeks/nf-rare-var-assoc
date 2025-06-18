from diagrams import Cluster, Diagram
from diagrams.c4 import Person, Container, Database, System, SystemBoundary, Relationship
from diagrams.aws.compute import ElasticContainerServiceService
from diagrams.aws.analytics import AmazonOpensearchService

graph_attr = {
    "splines": "spline",
}

# direction="LR",
with Diagram("Nextflow gene-level associacion pipeline - regenie", direction="LR", graph_attr=graph_attr):

    phenotype_file = ElasticContainerServiceService("Phenotype file")
    
    with Cluster("Reports and outputs"):
        eda_diagrams = AmazonOpensearchService("Input Data Characteristics Plots")
        pc_plot = AmazonOpensearchService("Principal Components Plots")
        html_report = AmazonOpensearchService("HTML Report with Mahattan and QQ plots")
        csv_reports = AmazonOpensearchService("CSV files with results")
      
    
    with Cluster("Data Preprocessing", direction="TB"):
        vcf_data = Container(
            name="Input data in VCF format",
            technology="",
            description="Before filtering for Regenie step 1",
        )
        
        plink_data = Container(
            name="Input data in Plink BED/BIM/FAM format",
            technology="",
            description="Before filtering for Regenie step 1",
        )
        
        phenotype_file >> Relationship("") >> plink_data
        vcf_data >> Relationship("") >> plink_data
        phenotype_file >> Relationship("") >> eda_diagrams
        vcf_data >> Relationship("") >> eda_diagrams
    
    with Cluster("Data Filtering for Regenie Step 1"):
        plink_f_filtering = Container(
            name="Plink",
            technology="",
            description="Inbreeding Coefficient Filtering",
        )
        
        plink_pca = Container(
            name="Plink",
            technology="",
            description="Principal Component Analysis",
        )
        plink_f_filtering >> Relationship("") >> plink_pca
        plink_pca >> Relationship("") >> pc_plot
            
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
        
        vcf_data >> Relationship("") >> rscript_prepare_annotations
        plink_data >> Relationship("") >> plink_pgen
        
        rscript_prepare_annotations >> Relationship("") >> regenie_step_2
        plink_pgen >> Relationship("") >> regenie_step_2
        
    build_reports = Container(
        name="RScript",
        technology="",
        description="Build reports",
    )
    
    
    plink_data >> Relationship("") >> plink_f_filtering
    
    regenie_step_2 >> Relationship("") >> build_reports
    eda_diagrams >> Relationship("") >> build_reports
    pc_plot >> Relationship("") >> build_reports
    
    build_reports >> Relationship("") >> html_report
    build_reports >> Relationship("") >> csv_reports
    
    
