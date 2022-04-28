import * as d3 from "d3";
import * as cons from "./constants";

export const zoomInN = () => {
    d3.select('svg')
      .transition()
      .call(zoomM.scaleBy, 2);
}
  
export const zoomOutN = () => {
    d3.select('svg')
      .transition()
      .call(zoomM.scaleBy, 0.5);
}

export const zoomM = 
d3.zoom()
.scaleExtent([0.5, 16])
// .translateExtent([])
.on('zoom', handleZoom);

function handleZoom(e) {
  d3.select('svg g')
  .attr('transform', e.transform);
}

export const resetZoom = () => {
    d3.select('svg')
       .transition()
       .call(zoomM.scaleTo, 1);
   }

