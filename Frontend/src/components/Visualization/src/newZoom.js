import * as d3 from "d3";
import * as cons from "./constants";

// export const zoomInN = () => {
//     d3.select('svg')
//       .transition()
//       .call(zoomM.scaleBy, 2);
// }
  
// export const zoomOutN = () => {
//     d3.select('svg')
//       .transition()
//       .call(zoomM.scaleBy, 0.5);
// }

// let zoomM = 
// d3.zoom()
// .scaleExtent([1, 16])
// .translateExtent([[0, 0], [cons.height, cons.width]])
// .on('zoom', handleZoom);

// function handleZoom(e) {
//   select('svg g')
//   .attr('transform', e.transform);
// }

// export const resetZoom = () => {
//     select('svg')
//        .transition()
//        .call(zoomM.scaleTo, 1);
//    }

let zoom = d3.zoom()
// .scaleExtent([1, 16])
// .translateExtent([[0, 0], [500, 500]])
.on("zoom", function () {
svg.attr("transform", d3.event.transform)
});

export const panZoom = () => {
  d3.select("svg")
  .call(zoom
    )
}

export const zoomInN = () => {
  zoom.scaleBy("svg".transition(), 2);
}

export const zoomOutN = () => {
  zoom.scaleBy("svg".transition(), 0.5);
}