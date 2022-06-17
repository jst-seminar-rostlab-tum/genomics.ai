import { increasedPointRadius, initialPointRadius } from './constants';

export const fadeIn = (selection, duration = 1000, delay = 0) => selection.attr('opacity', 0)
  .transition()
  .duration(duration)
  .delay(delay)
  .attr('opacity', 1);

export const fadeOut = (selection, duration = 1000, delay = 0) => selection.attr('opacity', 1)
  .transition()
  .duration(duration)
  .delay(delay)
  .attr('opacity', 0);

export const delay = (t) => new Promise((resolve) => setTimeout(resolve, t));

export const euclideanDist = (p, q) => Math.sqrt((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2);

export const affineComb = (l) => (1 - l) * initialPointRadius + l * increasedPointRadius;
