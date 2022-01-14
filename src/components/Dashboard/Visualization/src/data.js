import { scaleLinear, scaleOrdinal, min, max, mean, sum, quadtree } from 'd3';
import { width, height, colors, searchRadius } from './constants';
import { Cluster } from './cluster';

const prepare = ([metaData, umapData]) => {

  const clusterMapping = new Map();

  metaData
    .filter(({ CELLTYPE }) => isNaN(parseFloat(CELLTYPE)))
    .forEach(m => clusterMapping.set(parseInt(m.CLUSTER), m.CELLTYPE));

  const sampleData = umapData
    .map(d => ({
      id: d.cellType,
      x: +d.x,
      y: +d.y,
      cluster: parseInt(d.cluster)
    }))
    .filter(d => clusterMapping.has(d.cluster));

  return { clusterMapping, sampleData };
};

const getQuadTree = sampleData => {

  const tree = quadtree()
    .x(p => p.x)
    .y(p => p.y);

  tree.addAll(sampleData);
  return tree;
};

const getClusters = (clusterMapping, sampleData) => {

  const clusters = new Map();

  for (const p of sampleData) {
    const clusterName = clusterMapping.get(p.cluster);
    if (!clusters.has(clusterName))
      clusters.set(clusterName, []);
    clusters.get(clusterName).push([p.x, p.y]);
  }

  const names = clusters.keys();
  const colorScale = scaleOrdinal()
    .domain(...names)
    .range(colors);

  for (const name of clusters.keys()) {
    const points = Array.from(clusters.get(name));
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

    const { clusterMapping, sampleData } = prepare(data);
    this.tree = getQuadTree(sampleData);
    this.collection = getClusters(clusterMapping, sampleData);
    this.clusters = Array.from(this.collection.values());

    const { xMin, xMax, yMin, yMax } = getBounds(this.clusters);
    this.initialXDomain = scaleLinear()
        .domain([xMin, xMax])
        .range([0, width]);
    this.initialYDomain = scaleLinear()
        .domain([yMin, yMax])
        .range([height, 0]);
    this.xDomain = this.initialXDomain;
    this.yDomain = this.initialYDomain;

    this.totalCount = sum(this.clusters.map(c => c.count));

    const dists = this.clusters.map(c => c.dev);
    this.minDist = min(dists);
    this.maxDist = max(dists);
    this.avgDist = mean(dists);
  };

  getAll() {
    return this.clusters;
  };

  get(name) {
    return this.collection.get(name);
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