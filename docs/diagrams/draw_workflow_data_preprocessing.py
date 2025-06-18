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

def _format_description(description, width=30):
    """
    Formats the description string so it fits into the C4 nodes.

    It line-breaks the description so it fits onto exactly three lines. If there are more
    than three lines, all further lines are discarded and "..." inserted on the last line to
    indicate that it was shortened. This will also html-escape the description so it can
    safely be included in a HTML label.
    """
    wrapper = textwrap.TextWrapper(width=width, max_lines=4)
    lines = [html.escape(line) for line in wrapper.wrap(description)]
    # fill up with empty lines so it is always three
    lines += [""] * (3 - len(lines))
    return "<br/>".join(lines)

def _format_node_label(name, description, width=30):
    """Create a graphviz label string for a C4 node"""
    title = f'<font point-size="24"><b>{html.escape(name)}</b></font><br/>'
    text = f'<br/><font point-size="22">{_format_description(description, width)}</font>' if description else ""
    return f"<{title}{text}>"

def _format_edge_label(description):
    """Create a graphviz label string for a C4 edge"""
    wrapper = textwrap.TextWrapper(width=28, max_lines=3)
    lines = [html.escape(line) for line in wrapper.wrap(description)]
    text = "<br/>".join(lines)
    return f'<<font point-size="22">{text}</font>>'

def file_node(label: str):
	return ElasticContainerServiceService(f"\n\n\n\n{label}", fontsize="22", width = "1.3", height = "2.2", imagepos = "tc", labelloc="b", fixedsize="shape")

def plot_node(label: str):
	return AmazonOpensearchService(f"\n\n\n\n{label}", fontsize="22", width = "1.3", height = "2.2", imagepos = "tc", labelloc="b", fixedsize="shape")

def cluster(label: str):
	return Cluster(label, graph_attr={"height": "6.1", "fontsize": "26.0", "fixedsize": "false", "size": "5", "pad": "5.0"})

def container(name: str, description: str, width="4.2", height="1.7"):
	return Container(
		name=name,
		description=description,
		label = _format_node_label(name, description, width=float(width)*7),
		width = width,
		height = height
	)

# direction="LR",
with Diagram("Nextflow gene-level associacion pipeline - data preprocessing", direction="TB", graph_attr=graph_attr):
    
    phenotype_file = file_node("Phenotype file")
    final_pvcf_file = file_node("Final pVCF file")
    eda_diagrams = plot_node("Input Data Characteristics Plots")
    pc_plot = plot_node("Principal Components Plots")
    
    with cluster("Data Preprocessing"):
        
        exploratory_data_analysis = container(
            name="Python script",
            description="Exploratory Data Analysis",
        )
        
        bcftools_quality_filtering = container(
            name="BCFTools",
            description="1) Genotype quality and depth of coverage filtering          2) variant normalization, duplicates removal",
            height="1.9",
            width="4.5"
        )
        
        vep = container(
            name="VEP",
            description="Variant annotation",
        )
        
        extract_gnomad_allele_freq = container(
            name="BCFTools",
            description="Extract gnomad allele frequency data",
        )
        
        plink_missingness_per_pheno = container(
            name="Plink",
            description="Impute sex, missingness filtering per phenotype",
        )
        
        phenotype_file >> Relationship("") >> exploratory_data_analysis
        final_pvcf_file >> Relationship("") >> exploratory_data_analysis
        exploratory_data_analysis >> Relationship("") >> eda_diagrams
        
        final_pvcf_file >> Relationship("") >> bcftools_quality_filtering
        
        bcftools_quality_filtering >> Relationship("") >> vep
        
        vep >> Relationship("") >> extract_gnomad_allele_freq
        phenotype_file >> Relationship("") >> plink_missingness_per_pheno
        vep >> Relationship("") >> plink_missingness_per_pheno

    with cluster("Data Filtering for Regenie Step 1"):
        
        plink_3 = container(
            name="Plink",
            description="1) Minor allele frequency and count filtering                                  2) Hardy-Weinberg equilibrium deviations filtering",
            height="2.2",
            width="11.4"
        )
        
        with cluster("Inbreeding Coefficient Filtering"):
            plink_f_filtering = container(
                name="Plink",
                description="1) --indep-pairwise linkage equilibrium calculation                                2) --het F-coefficient calculation   3) Calculate inbreeding outliers in a Python script and then remove them using Plink",
                height="2.2",
                width="11.4"
            )
        
        with cluster("Principal Component Analysis"):
            
            plink_pca = container(
                name="Plink and Python script",
                description="1) Remove high-LD regions   2) Remove related individuals (--king-cutoff)                       3) Perform PCA   4) Draw Principal Components plot",
                height="2.2",
                width="11.4"
            )
            
            extract_gnomad_allele_freq >> Relationship("") >> plink_pca
            plink_pca >> Relationship("") >> pc_plot
        
        
        plink_missingness_per_pheno >> Relationship("") >> plink_3
        plink_3 >> Relationship("") >> plink_f_filtering
        plink_f_filtering >> Relationship("") >> plink_pca
        
    
    
