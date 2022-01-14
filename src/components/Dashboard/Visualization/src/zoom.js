import { width, height } from './constants';
import { easeCircle, zoom } from 'd3';

const transitionTime = 1000;
const maxZoomFactor = 16;

export const addZoom = (svg, handleZoom) => {

  const z = zoom()
    .scaleExtent([1, maxZoomFactor])
    .translateExtent([
        [0, 0],
        [width, height],
    ])
    .on('zoom', handleZoom);

  svg.call(z)
    .on('wheel.zoom', null);

  const setZoom = zoomFactor => {
    svg.transition()
      .ease(easeCircle)
      .duration(transitionTime)
      .call(z.scaleBy, zoomFactor);
  }

  return setZoom;
}
