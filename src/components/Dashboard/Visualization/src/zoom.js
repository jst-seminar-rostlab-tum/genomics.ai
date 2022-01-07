import * as d3 from 'd3';
import { width, height } from './constants';

const transitionTime = 1000;
const maxZoomFactor = 20;

const addZoom = (svg, handleZoom) => {
  const zoom = d3
    .zoom()
    .scaleExtent([1, maxZoomFactor])
    .translateExtent([
      [0, 0],
      [width, height],
    ])
    .on('zoom', handleZoom);

  svg.call(zoom);

  const setZoom = (zoomFactor) => {
    svg.transition()
      .ease(d3.easeCircle)
      .duration(transitionTime)
      .call(zoom.scaleBy, zoomFactor);
  };

  return setZoom;
};
export default addZoom;
