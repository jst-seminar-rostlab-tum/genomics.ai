import { addSmallCloseButton } from './buttons';
import { fadeIn, fadeOut } from './helpers';
import { addBarPlot } from './barPlot';
import { addPiePlot } from './piePlot';
import {
  width,
  height,
  padding,
  fill,
  outerPieRadius,
  labelHeight,
  barWidth,
  textHeight,
  textLabelHeight,
  barHeight,
  cornerRadius,
  offset
} from './constants';

const middle = labelHeight / 2;
const extendedWidth = 4 * padding + barWidth + 2 * outerPieRadius;
const extendedHeight = textHeight + textLabelHeight + barHeight + 2.5 * padding;

export class LabelFactory {

  constructor(container, collection) {

    this.container = container;
    this.collection = collection;
    this.extendedShown = false;
    this.allLabelsShown = false;
    this.selectedCluster = null;
    this.visibility = new Map();

    for (const c of collection.clusters)
      this.visibility.set(c.id, false);
  }

  showLabel(id, [px, py], text, duration = 1000) {

    const textLength = 12 * text.length;
    const labelWidth = textLength + padding;
    const x = this.collection.xDomain(px);
    const y = this.collection.yDomain(py) - middle;
    const flip = (width - x - extendedWidth) < 0;
    const diff = flip ? (1 - labelWidth - offset) : (offset - 1);

    const label = this.container.append('g')
      .attr('id', `label-${id}`)
      .attr('transform', `translate(${x},${y})`);

    label.append('polygon')
      .attr('points', [[0, middle], [offset, middle + 5], [offset, middle - 5]])
      .attr('transform', `translate(0,${flip ? 20 : 0}) rotate(${flip ? 180 : 0})`)
      .attr('fill', fill);

    const container = label.append('g')
      .attr('id', `label-${id}-container`)
      .attr('transform', `translate(${diff},0)`);

    container.append('rect')
      .attr('x', 0)
      .attr('y', 0)
      .attr('width', labelWidth)
      .attr('height', labelHeight)
      .attr('rx', cornerRadius)
      .attr('fill', fill);

    container.append('text')
      .attr('x', offset)
      .attr('y', middle + 5)
      .attr('stroke', 'lightgrey')
      .attr('fill', 'white')
      .attr('textLength', textLength)
      .attr('font-family', 'ui-sans-serif')
      .attr('font-weight', 600)
      .attr('font-style', 'italic')
      .text(text);

    fadeIn(label, duration);
  };

  removeLabel(id, duration = 1000) {
    const label = this.container.selectAll(`#label-${id}`);
    fadeOut(label, duration).remove();
  };

  showClusterLabel(cluster, duration = 600) {

    if (this.visibility.get(cluster.id))
      return;

    const id = cluster.id;
    const p = cluster.mean;
    const text = cluster.name;

    this.showLabel(id, p, text, duration);
    this.visibility.set(cluster.id, true);
  };

  removeClusterLabel(cluster, duration = 1000) {

    const shown = this.visibility.get(cluster.id);

    if (!shown || this.allLabelsShown)
      return;

    this.removeLabel(cluster.id, duration);
    this.selectedCluster = null;
    this.extendedShown = false;
    this.visibility.set(cluster.id, false);
  };

  showLabels(duration = 1000) {
    this.allLabelsShown = true;
    for (const c of this.collection.clusters)
      if (!this.visibility.get(c.id))
        this.showClusterLabel(c, duration);
  };

  removeLabels(duration = 1000) {

    this.allLabelsShown = false;
    this.extendedShown = false;
    this.selectedCluster = null;

    for (const c of this.collection.clusters)
      if (this.visibility.get(c.id))
        this.removeClusterLabel(c, duration);
  };

  retractClusterLabel(cluster) {

    const shown = this.visibility.get(cluster.id);

    if (!shown || !this.allLabelsShown)
      return;

    const text = cluster.name;
    const textLength = 12 * text.length;
    const labelWidth = textLength + 20;

    const px = this.collection.scalePoint(cluster.mean)[0];
    const shiftX = (width - px) <= extendedWidth;
    const diffX = shiftX ? (1 - offset - labelWidth) : (offset - 1);
    const labelContainer = this.container.select(`#label-${cluster.id}-container`);

    labelContainer.transition()
      .duration(800)
      .delay(500)
        .attr('transform', `translate(${diffX},0)`)

    labelContainer.select('rect')
      .transition()
        .duration(800)
        .delay(500)
          .attr('x', 0)
          .attr('y', 0)
          .attr('width', labelWidth)
          .attr('height', labelHeight)
          .attr('rx', cornerRadius)
          .attr('fill', fill);

    const content = labelContainer.selectAll('g');
    fadeOut(content, 500).remove();
    this.selectedCluster = null;
    this.extendedShown = false;
  };

  showDetailTag(cluster) {

    if (cluster === this.selectedCluster)
      return;

    if (this.extendedShown) {
      this.extendedShown = false;
      if (this.allLabelsShown)
        this.retractClusterLabel(this.selectedCluster);
      else
        this.removeClusterLabel(this.selectedCluster);
      this.showClusterLabel(cluster);
    }

    this.selectedCluster = cluster;
    this.extendedShown = true;

    const [px, py] = this.collection.scalePoint(cluster.mean);
    const shiftX = (width - px) <= extendedWidth;
    const shiftY = (height - py) <= extendedHeight;
    const diffX = shiftX ? (1 - offset - extendedWidth) : (offset - 1);
    const diffY = shiftY ? (height - py - extendedHeight) : 0;
    const labelContainer = this.container.select(`#label-${cluster.id}-container`);

    labelContainer.transition()
      .duration(500)
        .attr('transform', `translate(${diffX},${diffY})`);

    labelContainer.select('rect')
      .transition()
        .duration(500)
          .attr('width', extendedWidth)
          .attr('height', extendedHeight)

    addBarPlot(labelContainer, cluster, this.collection);
    addPiePlot(labelContainer, cluster, this.collection);

    const gClose = labelContainer.append('g')
      .attr('transform', `translate(${extendedWidth - 15 - 3},3)`)

    addSmallCloseButton(gClose);
    fadeIn(gClose, 0, 500);

    gClose.on('click', () => {
      this.selectedCluster = null;
      this.extendedShown = false;
      if (this.allLabelsShown)
        this.retractClusterLabel(cluster);
      else this.removeClusterLabel(cluster);
    });
  };

  allShown() {
    return this.allLabelsShown;
  };

  extended() {
    return this.extendedShown;
  };

  selected() {
    return this.selectedCluster;
  };
};