import * as d3 from "d3";

export const zoomInN = () => {

  d3.select('svg')
    .transition()
    .call(zoomM2.scaleBy, 1.75);
}

export const zoomOutN = () => {
  d3.select('svg')
    .transition()
    .call(zoomM2.scaleBy, 0.5);
}


export const zoomM2 =
  d3.zoom()
    .scaleExtent([0.5, 16])
    .on('zoom', handleZoom2);

function handleZoom2(e) {
  d3.select('svg g')
    .style("transform-origin", "50% 50% 0")
    .attr('transform', e.transform)

}

export const zoomM =
  d3.zoom()
    .scaleExtent([0.5, 16])
    .on('zoom', handleZoom);

function handleZoom(e) {
  d3.select('svg g')
    .style('transform-origin', null)
    .attr('transform', e.transform);
}

export const resetZoom = () => {
  d3.select('svg')
    .transition()
    .call(zoomM.scaleTo, 1);
}

