import * as d3 from "d3";


export let zoom = d3.zoom()
        .scaleExtent([0.25, 10])
        .on('zoom', handleZoom);

export function initZoom() {
        d3.select('svg g')
            .call(zoom);
    }

function handleZoom(e) {
    d3.select('svg g')
        .attr('transform', e.transform);
}

export function zoomIn() {
    d3.select('svg g')
        .transition()
        .call(zoom.scaleBy, 1.75);
}

export function zoomOut() {
    d3.select('svg g')
        .transition()
        .call(zoom.scaleBy, 0.5);
}

export function resetZoom() {
    d3.select('svg g')
        .transition()
        .call(zoom.scaleTo, 1);
}


/**export const zoomInN = () => {

  d3.select('svg g')
  .style("transform-origin", "20% 20% 0")
    .transition()
    .call(zoomM.scaleBy, 1.75);
}

export const zoomOutN = () => {
  d3.select('svg g')
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
  d3.select('svg g')
    .transition()
    .call(zoomM.scaleTo, 1);
}
*/
