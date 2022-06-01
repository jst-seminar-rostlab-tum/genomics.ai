import * as d3 from "d3";

// Zooms in a certain amount. Used for the zoom buttons
export const zoomInN = () => {

  d3.select('#vis_svg')
    .style("transform-origin", "20% 20% 0")
    .transition()
    .call(zoomM.scaleBy, 1.75);
}

// Zooms out a certain amount. Used for the zoom buttons
export const zoomOutN = () => {
  d3.select('#vis_svg')
    .style("transform-origin", "20% 20% 0")
    .transition()
    .call(zoomM.scaleBy, 0.5);
}

// Constructs the event-handler
export let zoomM =
  d3.zoom()
    .scaleExtent([0.5, 16])
    .on('zoom', handleZoom);

// Handles the zoom with the mouse wheel
export function handleZoom(e) {
  d3.select('svg g')
    .style("transform-origin", null)
    .attr('transform', e.transform);
}

// Resets the scale
export const resetZoom = () => {
  d3.select('#vis_svg')
    .transition()
    .call(zoomM.transform, d3.zoomIdentity);
}

