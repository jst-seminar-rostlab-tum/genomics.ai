import React from 'react';
import * as d3 from 'd3';
import umapData from './data/retinal_bipolar_viz/umap_data.txt';
import metaData from './data/retinal_bipolar_viz/meta_data.txt';
import { prepareContainer, addGroups } from './container';
import { width, height, margin } from './constants';
import { getDomains, prepareData } from './data';
import { addAxes, addAxisLabel } from './axis';
import addTooltip from './tooltip';
import addZoom from './zoom';
import styles from './Visualization.css';
import addLegend from './legend';

class Visualization extends React.Component {
  constructor(props) {
    super(props);
    this.myRef = React.createRef();
  }

  componentDidMount() {
    const zoomStep = 5;
    const initialRadius = 2;
    const labelPadding = 4;
    const middleX = margin.left + width / 2;
    const middleY = -margin.top - height / 2;
    const xPadding = margin.top - labelPadding;
    const yPadding = margin.left - labelPadding;
    const xAxisLabel = 'UMAP-1';
    const yAxisLabel = 'UMAP-2';

    const container = d3.select(this.myRef.current);
    let zoomFactor = 1;

    Promise.all([d3.csv(metaData), d3.csv(umapData)]).then((data) => {
      const { clusterMapping, pointsData } = prepareData(data);
      const { domX, domY, colorScale } = getDomains(clusterMapping, pointsData);
      const {
        svgDots, svgLegend, buttons, tooltip,
      } = prepareContainer(container);
      const { gX, gY, gDots } = addGroups(svgDots);
      const { xAxis, yAxis } = addAxes(gX, gY, domX, domY);

      addAxisLabel(svgDots, xAxisLabel, middleX, xPadding);
      addAxisLabel(svgDots, yAxisLabel, middleY, yPadding, 90);
      addLegend(svgLegend, clusterMapping, colorScale);

      const dots = gDots
        .selectAll('circle')
        .data(pointsData)
        .join('circle')
        .attr('cx', (d) => domX(d.x))
        .attr('cy', (d) => domY(d.y))
        .attr('fill', (d) => colorScale(d.cluster))
        .attr('r', initialRadius);

      addTooltip(tooltip, dots, clusterMapping);

      const handleZoom = ({ transform }) => {
        zoomFactor = transform.k;
        const dx = transform.rescaleX(domX);
        const dy = transform.rescaleY(domY);
        dots.attr('transform', transform);
        if (transform.k !== 1) dots.attr('r', initialRadius / (1 + Math.log2(zoomFactor)));
        gX.call(xAxis.scale(dx));
        gY.call(yAxis.scale(dy));
      };

      const setZoom = addZoom(svgDots, handleZoom);

      buttons.append('button')
        .text('ZOOM IN')
        .on('click', () => setZoom(zoomStep));

      buttons.append('button')
        .text('Zoom ouT'.toUpperCase())
        .on('click', () => setZoom(1 / zoomStep));

      buttons.append('button')
        .text('Reset zoom'.toUpperCase())
        .on('click', () => setZoom(1 / zoomFactor));
    });
  }

  render() {
    return <div ref={this.myRef} className={styles.container} />;
  }
}

export default Visualization;
