process GENERATE_TRACKING_REPORT {

    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/psuszynski/python_tools:1.0.0'

    input:
    tuple val(meta), path(tracking_files)
    
    output:
    tuple val(meta), path("*_sankey_report.html"), emit: report_html_file
    tuple val(meta), path("*_pipeline_report.txt"), emit: report_txt_file
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
import json
import sys
import os
import pandas as pd


# Read tracking files
tracking_files = "${tracking_files}".split(" ")
print(f"tracking_files = {tracking_files}")
#trace_file = sys.argv[-1]  # Last argument is trace.txt

# Parse trace.txt for execution metrics
metrics = {}
#if os.path.exists(trace_file):
#    try:
#        trace = pd.read_csv(trace_file, sep="\\t")
#        metrics = trace.groupby("name")[["realtime", "peak_rss"]].mean().to_dict()
#    except Exception as e:
#        print(f"Warning: Failed to parse {trace_file}: {e}")

# Initialize report
report_lines = [
    "Pipeline Report",
    "==============",
    ""
]

# Process each tracking file
for file in tracking_files:
    with open(file) as f:
        try:
            data = json.load(f)
        except ValueError as e:
            print(f"Error while reading file {file}, error: {e}")
            data = None
    
    if data is not None:
        # Extract and shorten process name
        full_name = data["process_name"]
        workflow_name = data["workflow_name"]
        short_name = full_name.split(":")[-1]
        
        # Extract variants and samples
        variants_in = data["inputs"]["variants"]
        samples_in = data["inputs"]["samples"]
        if "outputs" in data:
            variants_out = data["outputs"]["variants"]
            samples_out = data["outputs"]["samples"]
        else:
            variants_out = -1
            samples_out = -1
        
        # Extract parameters
        #params = ", ".join(f"{k}={v}" for k, v in data["parameters"].items())
        params = data["parameters"]

        # Extract and shorten predecessor
        predecessor = []
        if data["predecessor"] != "none":
            for pred_elem in data["predecessor"].split(" "):
                predecessor.append(pred_elem.split(":")[-1])
        else:
            predecessor = ["none"]
        predecessor = " ".join(predecessor)
        
        # Get execution metrics
        duration = metrics.get("realtime", {}).get(full_name, 0) / 1000  # ms to s
        memory = metrics.get("peak_rss", {}).get(full_name, 0) / 1024**2  # bytes to MB
        duration_str = f"{duration:.2f}" if duration > 0 else "N/A"
        memory_str = f"{memory:.2f}" if memory > 0 else "N/A"

        # Add to report
        report_lines.extend([
            f"Process: {short_name}",
            f"Workflow: {workflow_name}",
            f"Input Variants: {variants_in}",
            f"Output Variants: {variants_out}",
            f"Input Samples: {samples_in}",
            f"Output Samples: {samples_out}",
            f"Parameters: {params}",
            f"Predecessor: {predecessor}",
            f"Execution Time: {duration_str} s",
            f"Memory Usage: {memory_str} MB",
            "---------------",
            ""
        ])

# Write report
with open("${prefix}_pipeline_report.txt", "w") as f:
    f.write("\\n".join(report_lines))



import json
import os
from collections import defaultdict
import networkx as nx

# Function to load JSON files from a directory
def load_json_files():
    files = "${tracking_files}".split(" ")
    data = []
    for file in files:
        with open(file, 'r') as f:
            try:
                data.append(json.load(f))
            except ValueError as e:
                print(f"Error while reading file {file}, error: {e}")

    return data

# Function to build DAG and identify branches
def build_dag_and_branches(data):
    G = nx.DiGraph()
    process_dict = {item['process_name']: item for item in data}
    
    # Add nodes and edges
    for item in data:
        process_name = item['process_name']
        G.add_node(process_name, data=item)
        predecessors = item.get('predecessor', '').split()
        for pred in predecessors:
            if pred and pred in process_dict:
                G.add_edge(pred, process_name)
    
    # Identify disconnected subgraphs (branches)
    branches = list(nx.weakly_connected_components(G))
    return G, branches

# Function to prepare Sankey data for a branch
def prepare_sankey_data(G, branch, value_key):
    nodes = []
    links = []
    node_indices = {node: idx for idx, node in enumerate(branch)}
    
    # Add nodes
    for node in branch:
        if 'outputs' in G.nodes[node]['data'] and G.nodes[node]['data']['outputs'][value_key] >= 0:
            value = G.nodes[node]['data']['outputs'][value_key]
        else:
            value = G.nodes[node]['data']['inputs'][value_key]
        nodes.append({"name": node, "value": value})
    
    # Add links
    for src in branch:
        for dst in G.successors(src):
            if dst in branch:
                if 'outputs' in G.nodes[dst]['data'] and G.nodes[dst]['data']['outputs'][value_key] >= 0:
                    value = G.nodes[dst]['data']['outputs'][value_key]
                else:
                    value = G.nodes[dst]['data']['inputs'][value_key]
                links.append({
                    "source": node_indices[src],
                    "target": node_indices[dst],
                    "value": value
                })
    
    return {"nodes": nodes, "links": links}

# Function to generate HTML report
def generate_html_report(branches, G, output_file):
    html_content = '''
    <!DOCTYPE html>
    <html>
    <head>
        <title>Workflow Sankey Report</title>
        <script src="https://d3js.org/d3.v7.min.js"></script>
        <script src="https://unpkg.com/d3-sankey@0.12.3/dist/d3-sankey.min.js"></script>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            h1, h2 {{ color: #333; }}
            .sankey-plot {{ margin-bottom: 50px; }}
            .section {{ margin-top: 30px; }}
        </style>
    </head>
    <body>
        <h1>Workflow Analysis Report</h1>
        
        <script>
            function createSankey(data, containerId, width, height) {{
                const svg = d3.select(`#\${{containerId}}`)
                    .append("svg")
                    .attr("width", width)
                    .attr("height", height);
                
                const sankey = d3.sankey()
                    .nodeWidth(15)
                    .nodePadding(10)
                    .extent([[1, 1], [width - 1, height - 6]]);
                
                const {{nodes, links}} = sankey({{
                    nodes: data.nodes.map(d => ({{...d}})),
                    links: data.links.map(d => ({{...d}}))
                }});
                
                svg.append("g")
                    .selectAll("rect")
                    .data(nodes)
                    .enter()
                    .append("rect")
                    .attr("x", d => d.x0)
                    .attr("y", d => d.y0)
                    .attr("height", d => d.y1 - d.y0)
                    .attr("width", d => d.x1 - d.x0)
                    .attr("fill", "steelblue")
                    .append("title")
                    .text(d => `\${{d.name}}\\\\n\${{d.value}}`);
                
                svg.append("g")
                    .selectAll("path")
                    .data(links)
                    .enter()
                    .append("path")
                    .attr("d", d3.sankeyLinkHorizontal())
                    .attr("stroke", "gray")
                    .attr("stroke-width", d => Math.max(1, d.width))
                    .attr("fill", "none")
                    .style("opacity", 0.5)
                    .append("title")
                    .text(d => `\${{d.source.name}} → \${{d.target.name}}\\\\n\${{d.value}}`);
                
                svg.append("g")
                    .selectAll("text")
                    .data(nodes)
                    .enter()
                    .append("text")
                    .attr("x", d => d.x0 - 6)
                    .attr("y", d => (d.y1 + d.y0) / 2)
                    .attr("dy", "0.35em")
                    .attr("text-anchor", "end")
                    .text(d => d.name.split(":").pop())
                    .filter(d => d.x0 < width / 2)
                    .attr("x", d => d.x1 + 6)
                    .attr("text-anchor", "start");
            }}
        </script>

        <div class="section">
            <h2>Variants Flow</h2>
            {0}
        </div>
        
        <div class="section">
            <h2>Samples Flow</h2>
            {1}
        </div>
    </body>
    </html>
    '''
    
    # Generate plot divs for variants and samples
    variants_plots = ""
    samples_plots = ""
    
    for i, branch in enumerate(branches):
        # Variants Sankey data
        variants_data = prepare_sankey_data(G, branch, 'variants')
        variants_div_id = f"variants_sankey_{i}"
        variants_plots += f'<div class="sankey-plot"><h3>Branch {i+1}</h3><div id="{variants_div_id}"></div></div>'
        variants_plots += f'<script>createSankey({json.dumps(variants_data)}, "{variants_div_id}", 8000, 600);</script>'
        
        # Samples Sankey data
        samples_data = prepare_sankey_data(G, branch, 'samples')
        samples_div_id = f"samples_sankey_{i}"
        samples_plots += f'<div class="sankey-plot"><h3>Branch {i+1}</h3><div id="{samples_div_id}"></div></div>'
        samples_plots += f'<script>createSankey({json.dumps(samples_data)}, "{samples_div_id}", 8000, 600);</script>'
    
    # Write HTML file
    with open(output_file, 'w') as f:
        f.write(html_content.format(variants_plots, samples_plots))


# Load data
data = load_json_files()
output_file = "${prefix}_sankey_report.html"

# Build DAG and identify branches
G, branches = build_dag_and_branches(data)
print(f"len(branches) = {len(branches)}")

# Generate HTML report
generate_html_report(branches, G, output_file)
print(f"Report generated: {output_file}")


# Write versions.yml
with open('versions.yml', 'w') as f:
    f.write('${task.process}:\\n')
    f.write(f'    python: {sys.version.split()[0]}\\n')
    f.write(f'    d3js: v7\\n')
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sankey_report.html
    touch ${prefix}_pipeline_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS
    """
}
