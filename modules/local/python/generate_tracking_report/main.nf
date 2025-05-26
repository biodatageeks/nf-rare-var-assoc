process GENERATE_TRACKING_REPORT {

    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seaborn:0.13.2':
        'biocontainers/seaborn:0.13.2' }"

    input:
    tuple val(meta), path(tracking_files)
    
    output:
    tuple val(meta), path("pipeline_report.html"), emit: report_html_file
    tuple val(meta), path("pipeline_report.txt"), emit: report_txt_file
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python
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
        predecessor = data["predecessor"].split(":")[-1] if data["predecessor"] != "none" else "none"
        # Get execution metrics
        duration = metrics.get("realtime", {}).get(full_name, 0) / 1000  # ms to s
        memory = metrics.get("peak_rss", {}).get(full_name, 0) / 1024**2  # bytes to MB
        duration_str = f"{duration:.2f}" if duration > 0 else "N/A"
        memory_str = f"{memory:.2f}" if memory > 0 else "N/A"

        # Add to report
        report_lines.extend([
            f"Process: {short_name}",
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
with open("pipeline_report.txt", "w") as f:
    f.write("\\n".join(report_lines))



nodes = []
edges = []
for file in tracking_files:
    with open(file) as f:
        try:
            data = json.load(f)
        except ValueError as e:
            print(f"Error while reading file {file}, error: {e}")
            data = None
    
    if data is not None:
        # Create node label
        # params = "\\n".join(f"{k}={v}" for k, v in data["parameters"].items())
        params = data["parameters"]
        label = (f"{data['process_name']}\\n"
                f"Variants: {data['inputs']['variants']}→{data.get('outputs', {'variants': -1})['variants']}\\n"
                f"Samples: {data['inputs']['samples']}→{data.get('outputs', {'samples': -1})['samples']}\\n"
                f"Parameters:\\n{params}")
        nodes.append({"id": data["process_name"], "label": label})
        if data["predecessor"] != "none":
            edges.append({"source": data["predecessor"], "target": data["process_name"]})

# HTML template with D3.js
html_template = '''
<!DOCTYPE html>
<html>
<head>
    <title>Pipeline Flow Diagram</title>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <style>
        body { font-family: Arial, sans-serif; }
        .node rect {
            fill: #f0f8ff;
            stroke: #4682b4;
            stroke-width: 1.5px;
        }
        .node text {
            font-size: 12px;
            text-anchor: middle;
        }
        .edgePath path {
            stroke: #333;
            stroke-width: 1.5px;
            fill: none;
        }
        .tooltip {
            position: absolute;
            background: #333;
            color: white;
            padding: 5px;
            border-radius: 3px;
            pointer-events: none;
        }
        svg { border: 1px solid #ccc; background: #fff; }
    </style>
</head>
<body>
    <h1>Pipeline Flow Diagram</h1>
    <svg width="1200" height="1000"></svg>
    <div id="tooltip" class="tooltip" style="opacity: 0;"></div>
    <script>
        try {
            const nodes = %s;
            const edges = %s;

            console.log("Nodes:", nodes);
            console.log("Edges:", edges);

            const svg = d3.select("svg"),
                width = +svg.attr("width"),
                height = +svg.attr("height");

            // Function to measure text width
            function getTextWidth(text, fontSize = 12) {
                const tempSvg = d3.select("body").append("svg");
                const tempText = tempSvg.append("text")
                    .attr("font-size", fontSize)
                    .text(text);
                const width = tempText.node().getBBox().width;
                tempSvg.remove();
                return width;
            }

            // Calculate node dimensions
            nodes.forEach(d => {
                const lines = d.label.split('\\\\n');
                d.width = Math.max(...lines.map(line => getTextWidth(line, 12))) + 20; // Add padding
                // d.height = lines.length * 18 + 20; // 18px per line + padding
                d.height = lines.length * 24 + 100;
            });

            // Force simulation
            const simulation = d3.forceSimulation(nodes)
                .force("link", d3.forceLink(edges).id(d => d.id).distance(150))
                .force("charge", d3.forceManyBody().strength(-500))
                .force("center", d3.forceCenter(width / 2, height / 2))
                .force("x", d3.forceX().strength(0.1))
                .force("y", d3.forceY().strength(0.1));

            // Add links
            const link = svg.append("g")
                .attr("class", "links")
                .selectAll("line")
                .data(edges)
                .enter().append("line")
                .attr("class", "edgePath")
                .attr("stroke", "#333")
                .attr("stroke-width", 1.5)
                .attr("marker-end", "url(#arrow)");

            // Add nodes
            const node = svg.append("g")
                .attr("class", "nodes")
                .selectAll("g")
                .data(nodes)
                .enter().append("g")
                .attr("class", "node")
                .call(d3.drag()
                    .on("start", dragstarted)
                    .on("drag", dragged)
                    .on("end", dragended))
                .on("mouseover", showTooltip)
                .on("mouseout", hideTooltip);

            node.append("rect")
                .attr("width", d => d.width)
                .attr("height", d => d.height)
                .attr("x", d => -d.width / 2)
                .attr("y", d => -d.height / 2);

            node.append("text")
                .selectAll("tspan")
                .data(d => d.label.split('\\\\n'))
                .enter().append("tspan")
                .attr("x", 0)
                .attr("dy", (d, i) => i * 18)
                .text(d => d);

            // Arrowheads
            svg.append("defs").append("marker")
                .attr("id", "arrow")
                .attr("viewBox", "0 -5 10 10")
                .attr("refX", 15)  // Adjusted for better arrow placement
                .attr("refY", 0)
                .attr("markerWidth", 8)
                .attr("markerHeight", 8)
                .attr("orient", "auto")
                .append("path")
                .attr("d", "M0,-5L10,0L0,5")
                .attr("fill", "#333");

            // Update positions
            simulation.on("tick", () => {
                link
                    .attr("x1", d => d.source.x)
                    .attr("y1", d => d.source.y)
                    .attr("x2", d => d.target.x)
                    .attr("y2", d => d.target.y);
                node.attr("transform", d => {
                    d.x = Math.max(d.width / 2, Math.min(width - d.width / 2, d.x));
                    d.y = Math.max(d.height / 2, Math.min(height - d.height / 2, d.y));
                    return `translate(\${d.x},\${d.y})`;
                });
            });

            function dragstarted(event, d) {
                if (!event.active) simulation.alphaTarget(0.3).restart();
                d.fx = d.x;
                d.fy = d.y;
            }
            function dragged(event, d) {
                d.fx = event.x;
                d.fy = event.y;
            }
            function dragended(event, d) {
                if (!event.active) simulation.alphaTarget(0);
                d.fx = null;
                d.fy = null;
            }

            const tooltip = d3.select("#tooltip");
            function showTooltip(event, d) {
                tooltip.style("opacity", 1)
                    .html(`Full Name: \${d.full_name}`)
                    .style("left", (event.pageX + 10) + "px")
                    .style("top", (event.pageY - 10) + "px");
            }
            function hideTooltip() {
                tooltip.style("opacity", 0);
            }

        } catch (error) {
            console.error("D3.js Error:", error);
            alert("Error rendering diagram: " + error.message);
        }
    </script>
</body>
</html>
'''

# Write HTML report
with open("pipeline_report.html", "w") as f:
    f.write(html_template % (json.dumps(nodes), json.dumps(edges)))

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
    touch ${prefix}_${out_name_part}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS
    """
}