import { scaleLinear, scaleOrdinal, min, max, mean, sum, quadtree } from 'd3';
import { colors, searchRadius } from './constants';
import { Cluster } from './cluster';

const getQuadTree = data => {

  const tree = quadtree()
    .x(p => +p.x)
    .y(p => +p.y);

  tree.addAll(data);
  return tree;
};

const getClusters = data => {

  const clusters = new Map();
  data
    .filter(({ celltype }) => isNaN(parseFloat(celltype)))
    .forEach(p => {
      const clusterName = p.celltype;
      if (!clusters.has(clusterName))
        clusters.set(clusterName, []);
      clusters.get(clusterName).push([+p.x, +p.y]);
    });

  const names = clusters.keys();
  const colorScale = scaleOrdinal()
    .domain(...names)
    .range(colors);

  for (const name of clusters.keys()) {
    const points = clusters.get(name);
    const color = colorScale(name);
    clusters.set(name, new Cluster(points, name, color));
  }

  return clusters;
};

const getBounds = clusters => {

  const xMin = min(clusters.map(c => c.getMinX()));
  const xMax = max(clusters.map(c => c.getMaxX()));
  const yMin = min(clusters.map(c => c.getMinY()));
  const yMax = max(clusters.map(c => c.getMaxY()));

  return { xMin, xMax, yMin, yMax };
};

export class ClusterCollection {

  constructor(data) {

    this.tree = getQuadTree(data);
    this.collection = getClusters(data);
    this.clusters = Array.from(this.collection.values());

    const { xMin, xMax, yMin, yMax } = getBounds(this.clusters);
    this.xMin = xMin;
    this.xMax = xMax;
    this.yMin = yMin;
    this.yMax = yMax;

    const dists = this.clusters.map(c => c.dev);
    this.minDist = min(dists);
    this.maxDist = max(dists);
    this.avgDist = mean(dists);

    this.totalCount = sum(this.clusters.map(c => c.count));
  };

  getAll() {
    return this.clusters;
  };

  get(name) {
    return this.collection.get(name);
  };

  setRange([width, height]) {
    this.width = width;
    this.height = height;
    this.initialXDomain = scaleLinear()
        .domain([this.xMin, this.xMax])
        .range([0, width]);
    this.initialYDomain = scaleLinear()
        .domain([this.yMin, this.yMax])
        .range([height, 0]);
    this.xDomain = this.initialXDomain;
    this.yDomain = this.initialYDomain;
  };

  rescaleDomain(transform) {
    this.xDomain = transform.rescaleX(this.initialXDomain);
    this.yDomain = transform.rescaleY(this.initialYDomain);
  };

  findClosestSample([x, y]) {
    return this.tree.find(x, y, searchRadius);
  };

  scalePoint([x, y]) {
    const xScaled = this.xDomain(x);
    const yScaled = this.yDomain(y);
    return [xScaled, yScaled];
  };

  unscalePoint([xScaled, yScaled]) {
    const x = this.xDomain.invert(xScaled);
    const y = this.yDomain.invert(yScaled);
    return [x, y];
  };

  scalePoints(points) {
    return points.map(p => this.scalePoint(p));
  };
};