import { margin } from "./constants";

const addColorLabels = (selection, clusterMapping, colorScale) => {
  const clusterLabels = clusterMapping.keys();
  selection
    .selectAll("colorLabels")
    .data(clusterLabels)
    .join("circle")
    .attr("cx", 10)
    .attr("cy", (d, i) => margin.top + 6 + i * 20)
    .attr("r", 6)
    .style("fill", c => colorScale(c));
}

const addTextLabels = (selection, clusterMapping) => {
  const clusterLabels = clusterMapping.keys();
  selection
    .selectAll("textLabels")
    .data(clusterLabels)
    .join("text")
    .attr("x", 20)
    .attr("y", (d, i) => margin.top + 9 + i * 20)
    .style("font-size", "10px")
    .text(c => clusterMapping.get(c));
}

export const addLegend = (svgLegend, clusterMapping, colorScale) => {
  const gLegend = svgLegend.append("g");
  addColorLabels(gLegend, clusterMapping, colorScale);
  addTextLabels(gLegend, clusterMapping);
  return gLegend;
}