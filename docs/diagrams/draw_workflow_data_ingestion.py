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
    wrapper = textwrap.TextWrapper(width=30, max_lines=4)
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
    return ElasticContainerServiceService(f"\n\n\n\n{label}", fontsize="22", width = "1.3", height = "2.2", imagepos = "tc", labelloc="b", fixedsize="shape")

def plot_node(label: str):
    return AmazonOpensearchService(f"\n\n\n\n{label}", fontsize="22", width = "1.3", height = "2.2", imagepos = "tc", labelloc="b", fixedsize="shape")

def cluster(label: str):
    return Cluster(label, graph_attr={"height": "6.1", "fontsize": "26.0", "fixedsize": "false", "size": "5", "pad": "5,10"})

def container(name: str, description: str, width="4.2", height="1.7", bgcolor=None):
    if bgcolor is not None:
        return Container(
            name=name,
            description=description,
            label = _format_node_label(name, description),
            width = width,
            height = height,
            fillcolor = bgcolor
        )
    else:
        return Container(
            name=name,
            description=description,
            label = _format_node_label(name, description),
            width = width,
            height = height
        )
        


# direction="LR",
with Diagram("Nextflow gene-level associacion pipeline - data ingestion", direction="LR", graph_attr=graph_attr):

    phenotype_file = file_node("Phenotype file")
    
    with cluster("pVCF Data Ingestion"):
        fasta_file = file_node("FASTA reference genome")
        bam_files = file_node("BAM/CRAM files")
        gvcf_files = file_node("gVCF files")
        
        deep_variant = container(
            name="DeepVariant",
            description="",
        )

        merge_gvcfs = container(
            name="Merge gvcfs",
            description="",
        )
        
        with cluster("Alternative Joint Variant Calling"):
            gl_nexus = container(
                name="GLNexus",
                description="",
            )
            
            new_ml_joint_variant_calling = container(
                name="New ML Joint Variant calling",
                description="developed as part of this project",
                bgcolor="darkorchid3",
                width="5.0"
            )

        bam_files >> Relationship("") >> deep_variant
        fasta_file >> Relationship("") >> deep_variant
        
        gvcf_files >> Relationship("") >> merge_gvcfs
        deep_variant >> Relationship("") >> merge_gvcfs
        
        merge_gvcfs >> Relationship("", color = "firebrick4") >> gl_nexus
        merge_gvcfs >> Relationship("", color = "firebrick4") >> new_ml_joint_variant_calling
        
        with cluster("Alternative Input"):
            input_pvcf_file = file_node("Input pVCF file")
            final_pvcf_file = file_node("Final pVCF file")
            
            gl_nexus >> Relationship("", color = "firebrick4") >> final_pvcf_file
            new_ml_joint_variant_calling >> Relationship("", color = "firebrick4") >> final_pvcf_file
            input_pvcf_file >> Relationship("", color = "firebrick4") >> final_pvcf_file
           
    
