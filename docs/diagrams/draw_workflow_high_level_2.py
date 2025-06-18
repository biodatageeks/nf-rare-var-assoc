import html
import textwrap

from diagrams import Cluster, Diagram, Edge
from diagrams.c4 import Person, Container, Database, System, SystemBoundary, Relationship
from diagrams.aws.compute import ElasticContainerServiceService
from diagrams.aws.analytics import AmazonOpensearchService

graph_attr = {
    "fontsize": "24.0",
    "splines": "spline",
    "label": "",
    "pad": "0"
}

node_attr = {
    "fontsize": "24.0"
}

def _format_description(description):
    """
    Formats the description string so it fits into the C4 nodes.

    It line-breaks the description so it fits onto exactly three lines. If there are more
    than three lines, all further lines are discarded and "..." inserted on the last line to
    indicate that it was shortened. This will also html-escape the description so it can
    safely be included in a HTML label.
    """
    wrapper = textwrap.TextWrapper(width=24, max_lines=3)
    lines = [html.escape(line) for line in wrapper.wrap(description)]
    # fill up with empty lines so it is always three
    lines += [""] * (3 - len(lines))
    return "<br/>".join(lines)

def _format_node_label(name, description):
    """Create a graphviz label string for a C4 node"""
    title = f'<font point-size="24"><b>{html.escape(name)}</b></font><br/>'
    text = f'<br/><font point-size="22">{_format_description(description)}</font>' if description else ""
    return f"<{title}{text}>"

def _format_edge_label(description):
    """Create a graphviz label string for a C4 edge"""
    wrapper = textwrap.TextWrapper(width=28, max_lines=3)
    lines = [html.escape(line) for line in wrapper.wrap(description)]
    text = "<br/>".join(lines)
    return f'<<font point-size="22">{text}</font>>'

def file_node(label: str):
	#return ElasticContainerServiceService(f"\n\n\n\n{label}", fontsize="18", width = "2.2", height = "1.1", imagepos = "tc", labelloc="t")
    return ElasticContainerServiceService(f"\n\n\n\n{label}", fontsize="22", width = "1.3", height = "2.2", imagepos = "tc", labelloc="b", fixedsize="shape")

def plot_node(label: str):
	#return AmazonOpensearchService(f"\n\n\n\n{label}", fontsize="18", width = "2.2", height = "1.1", imagepos = "tc", labelloc="t")
	return AmazonOpensearchService(f"\n\n\n\n{label}", fontsize="22", width = "1.3", height = "2.2", imagepos = "tc", labelloc="b", fixedsize="shape")

def cluster(label: str):
	return Cluster(label, graph_attr={"height": "6.1", "fontsize": "26.0", "fixedsize": "false", "size": "5", "pad": "5,10"})

def container(name: str, description: str, width="4.2", height="1.7"):
	return Container(
		name=name,
		description=description,
		label = _format_node_label(name, description),
		width = width,
		height = height
	)

# direction="LR",
with Diagram("Nextflow gene-level associacion pipeline - high level", direction="TB", graph_attr=graph_attr, node_attr=node_attr):
    
    with cluster("Inputs"):
        phenotype_file = file_node("Phenotype file")
        fasta_file = file_node("FASTA reference")
        bam_files = file_node("BAM/CRAM files")
        gvcf_files = file_node("gVCF files")
        input_pvcf_file = file_node("Input pVCF file")
    
    with cluster("Reports and outputs"):
        eda_diagrams = plot_node("Input Data Characteristics\n Plots")
        pc_plot = plot_node("Principal Components Plots")
        html_report = plot_node("HTML Report with \nMahattan and QQ plots")
        csv_reports = plot_node("CSV files with results")  
    
    with cluster("Data Ingestion and Preprocessing"):
        vcf_data = container(
            name="Input data in VCF format",
            description="Before filtering for Regenie step 1"
        )
        
        plink_data = container(
            name="Input data in Plink format",
            description="Before filtering for Regenie step 1",
            width="4.4"
        )
        
        fasta_file >> Relationship("") >> vcf_data
        bam_files >> Relationship("") >> vcf_data
        gvcf_files >> Relationship("") >> vcf_data
        input_pvcf_file >> Relationship("") >> vcf_data
        
        phenotype_file >> Relationship("") >> plink_data
        #vcf_data >> Relationship("") >> plink_data
        phenotype_file >> Relationship("") >> eda_diagrams
        vcf_data >> Relationship("") >> eda_diagrams
    
    with cluster("Data Filtering for Regenie Step 1"):
        plink_f_filtering = container(
            name="Plink",
            description="Inbreeding Coefficient Filtering",
        )
        
        plink_pca = container(
            name="Plink",
            description="Principal Component Analysis",
        )
        #plink_f_filtering >> Relationship("") >> plink_pca
        plink_pca >> Relationship("") >> pc_plot
            
    with cluster("Regenie"):
        
        regenie_step_1 = container(
            name="Regenie",
            description="Step 1",
            width="3.0"
        )
        
        
        rscript_prepare_annotations = container(
            name="RScript",
            description="Assign masks",
            width="3.0"
        )
        
        plink_pgen = container(
            name="Plink",
            description="Make PGEN/PVAR/PSAM files for optimal Regenie performance",
            width="4.0"
        )
    
        regenie_step_2 = container(
            name="Regenie",
            description="Step 2",
            width="3.0"
        )
    
        phenotype_file >> Relationship("") >> regenie_step_1
        plink_f_filtering >> Relationship("") >> regenie_step_1
        plink_pca >> Edge(
            style = "dashed",
            color = "gray60",
            label = _format_edge_label("covariates file")
        ) >> regenie_step_1
        
        vcf_data >> Relationship("") >> rscript_prepare_annotations
        plink_data >> Relationship("") >> plink_pgen
        
        rscript_prepare_annotations >> Relationship("") >> regenie_step_2
        plink_pgen >> Relationship("") >> regenie_step_2
        regenie_step_1 >> Relationship("") >> regenie_step_2
        
    build_reports = container(
        name="RScript",
        description="Build reports",
        width="3.0"
    )
    
    plink_data >> Relationship("") >> plink_f_filtering
    
    regenie_step_2 >> Relationship("") >> build_reports
    eda_diagrams >> Relationship("") >> build_reports
    pc_plot >> Relationship("") >> build_reports
    
    build_reports >> Relationship("") >> html_report
    build_reports >> Relationship("") >> csv_reports
    
    
