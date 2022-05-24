import * as d3 from "d3";

export const zoomInN = () => {

  d3.select('#vis_svg')
  .style("transform-origin", "20% 20% 0")
    .transition()
    .call(zoomM.scaleBy, 1.75);
}

export const zoomOutN = () => {
  d3.select('#vis_svg')
  .style("transform-origin", "20% 20% 0")
    .transition()
    .call(zoomM.scaleBy, 0.5);
}

export let zoomM =
  d3.zoom()
    .scaleExtent([0.5, 16])
    .on('zoom', handleZoom);

export function handleZoom(e) {
  d3.select('svg g')
    .style("transform-origin", null)
    .attr('transform', e.transform);
}

export const resetZoom = () => {
  d3.select('#vis_svg')
    .transition()
    .call(zoomM.scaleTo, 1);
}

