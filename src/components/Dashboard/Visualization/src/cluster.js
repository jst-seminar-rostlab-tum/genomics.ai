import { mean, extent } from 'd3';
import { euclideanDist } from './helpers';

export class Cluster {
  constructor(points, name, color) {
    const xs = points.map((p) => p[0]);
    const ys = points.map((p) => p[1]);
    const m = [mean(xs), mean(ys)];
    const devs = points.map((p) => euclideanDist(p, m));

    this.name = name;
    this.id = name.replaceAll(/[/,_, ]/g, '').toLowerCase();
    this.color = color;
    this.count = points.length;
    this.points = points;
    this.mean = m;
    this.dev = mean(devs);
    this.xExtent = extent(points.map((p) => +p[0]));
    this.yExtent = extent(points.map((p) => +p[1]));
  }

  getMinX() {
    return this.xExtent[0];
  }

  getMaxX() {
    return this.xExtent[1];
  }

  getMinY() {
    return this.yExtent[0];
  }

  getMaxY() {
    return this.yExtent[1];
  }
}
