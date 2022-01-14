import { Delaunay } from 'd3';
import { ClusterCollection } from './data';
import { addControlButtons } from './buttons';
import { LabelFactory } from './labels';
import { addZoom } from './zoom';
import { affineComb, delay } from './helpers';
import {
  width,
  height,
  centroidRadius,
  centroidBorderWidth,
  cellBorderWidth,
  padding,
  controlX,
  initialPointRadius,
  fill,
  zoomStep
} from './constants';

const addCanvas = container => {
  return container.append('canvas')
    .attr('width', width)
    .attr('height', height)
    .node()
    .getContext('2d');
};

const addSvg = container => {
  return container.append('svg')
    .attr('width', width)
    .attr('height', height)
    .attr('transform', `translate(0,${-height - 4})`);
    // .attr('preserveAspectRatio', 'xMinYMin meet')
    // .attr('viewBox', '0 0 600 400')
    // .classed('svg-content-responsive', true);
};

const addGroup = (svg, id, x = 0, y = 0) => {
  return svg.append('g')
    .attr('id', `${id}`)
    .attr('transform', `translate(${x},${y})`);
};

export class UmapVisualization {

  constructor(container, data) {

    this.context = addCanvas(container);
    this.svg = addSvg(container);

    this.gCentroids = addGroup(this.svg, 'centroids');
    this.gCells = addGroup(this.svg, 'cells');
    this.gLabels = addGroup(this.svg, 'labels');
    this.gControls = addGroup(this.svg, 'controls', controlX, padding);

    this.collection = new ClusterCollection(data);
    this.labels = new LabelFactory(this.svg, this.collection);

    const addCentroids = () => {

      const centroids = this.gCentroids.selectAll('centroids')
        .data(this.collection.clusters)
        .join('circle')
          .attr('cx', c => this.collection.xDomain(c.mean[0]))
          .attr('cy', c => this.collection.yDomain(c.mean[1]))
          .attr('r', centroidRadius)
          .attr('stroke', c => c.color)
          .attr('stroke-width', centroidBorderWidth)
          .attr('fill', c => c.color)
          .attr('fill-opacity', 0.7);

      return centroids;
    };

    const addCells = () => {

      const means = this.collection.clusters.map(c => c.mean);
      const scaledMeans = this.collection.scalePoints(means);
      const voronoi = Delaunay.from(scaledMeans)
        .voronoi([0, 0, width, height]);
      const cells = this.gCells.selectAll('cell')
        .data(this.collection.clusters)
        .join('polygon')
          .attr('id', c => `cell-${c.id}`)
          .attr('points', (c, i) => voronoi.cellPolygon(i))
          .attr('fill', c => c.color)
          .attr('opacity', 0.1)
          .attr('stroke', 'black')
          .attr('stroke-width', cellBorderWidth);

      return cells;
    };

    this.drawPoints(initialPointRadius);
    this.centroids = addCentroids();
    this.cells = addCells();

    this.cells.on('mouseover', (e, d) => {
      if (!this.labels.extended())
        this.labels.showClusterLabel(d);
      this.highlightCell(d);
    });

    this.cells.on('mouseleave', (e, d) => {
      if (!this.labels.extended())
        this.labels.removeClusterLabel(d);
      this.resetCell(d);
    });

    this.cells.on('click', (e, d) => {
        this.labels.showDetailTag(d);
    });

    const {
      zoomIn,
      zoomOut,
      toggleLabels,
      reset
    } = addControlButtons(this.gControls);

    this.zoomFactor = 1;
    this.counter = 0;
    this.zoomIn = zoomIn;
    this.zoomOut = zoomOut;
    this.toggleLabels = toggleLabels;
    this.reset = reset;

    this.zoomIn.on('click', () => {
      this.setZoom(zoomStep);
    });

    this.zoomOut.on('click', () => {
      this.setZoom(1 / zoomStep);
    });

    this.toggleLabels.on('mouseover', () => {
      if (!this.labels.allShown())
        this.toggleLabels.select('rect')
          .transition()
            .attr('fill', 'dimgray');
    });

    this.toggleLabels.on('mouseleave', () => {
      if (!this.labels.allShown())
        this.toggleLabels.select('rect')
          .transition()
            .attr('fill', fill);
    });

    this.toggleLabels.on('click', () => {

      if (!this.labels.allShown()) {
        this.labels.showLabels();
        this.toggleLabels.select('rect')
          .transition()
            .attr('fill', 'grey');
        return;
      }

    this.labels.removeLabels();
      this.toggleLabels.select('rect')
        .transition()
          .attr('fill', 'dimgray');
    });

    this.reset.on('click', () => {
      this.setZoom(1 / this.zoomFactor);
    });

    this.currentClosest = null;

    this.svg.on('mousemove', e => {

      if (this.zoomFactor >= 10) {

        const p = this.collection.unscalePoint([e.clientX, e.clientY]);
        const closest = this.collection.findClosestSample(p);

        if (closest !== this.currentClosest) {

          if (closest !== undefined) {
            const id = closest.id.replaceAll('_', '-');
            const p = [closest.x, closest.y];
            this.labels.showLabel(id, p, id);
          }

          if (this.currentClosest) {
            const id = this.currentClosest.id.replaceAll('_', '-');
            this.labels.removeLabel(id);
          }

          this.currentClosest = closest;
        };
      }
    });

    const handleZoom = ({ transform }) => {

      if (transform.k === 1)
        return;

      const allShown = this.labels.allShown();
      this.collection.rescaleDomain(transform);
      this.labels.removeLabels(0);
      this.clearCanvas();
      this.cells.attr('transform', transform);
      this.centroids.attr('transform', transform);

      if (transform.k !== this.zoomFactor) {

        this.zoomFactor = transform.k;
        this.centroids
          .attr('r', centroidRadius / this.zoomFactor)
          .attr('stroke-width', centroidBorderWidth / this.zoomFactor);
        this.cells
          .attr('stroke-width', cellBorderWidth / this.zoomFactor);

        if (allShown)
          this.labels.showLabels();

      } else {

        if (allShown) {
          this.counter++;
          delay(100).then(() => {
            this.counter--;
            if (this.counter === 0)
              this.labels.showLabels();
          });
        }
      }

      const r = affineComb(this.zoomFactor / 10);
      this.drawPoints(r);
    };

    this.setZoom = addZoom(this.svg, handleZoom);
  };

  drawPoints(r) {
    for (const c of this.collection.clusters) {
      this.context.fillStyle = c.color;
      this.context.globalAlpha = 0.3;
      for (const p of c.points) {
        const [x, y] = this.collection.scalePoint(p);
        this.context.beginPath();
        this.context.arc(x, y, r, 0, r * Math.PI);
        this.context.fill();
      }
    }
  };

  clearCanvas() {
    this.context.clearRect(0, 0, width, height);
  };

  highlightCell(cluster) {
    this.svg.select(`#cell-${cluster.id}`)
      .transition()
        .duration(500)
          .attr('opacity', 0.2);
  };

  resetCell(cluster) {
    this.svg.select(`#cell-${cluster.id}`)
      .transition()
        .duration(1000)
          .attr('opacity', 0.1);
  };
};