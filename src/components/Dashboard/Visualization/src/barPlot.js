import { padding, textHeight, textLabelHeight, barWidth, barHeight } from './constants';
import { scaleBand, scaleLinear, format, axisBottom, axisLeft } from 'd3';

const tickSize = 20;
const tickPadding = 2;
const innerTickSize = 5;
const outerTickSize = 0;
const labelFormatter = format('.2');

export const addBarPlot = (labelContainer, cluster, clusterCollection) => {

  const groups = ['', 'min', 'Ø', 'max'];
  const minDist = clusterCollection.minDist;
  const maxDist = clusterCollection.maxDist;
  const avgDist = clusterCollection.avgDist;
  const barData = [cluster.dev, minDist, avgDist, maxDist];

  const barPlot = labelContainer.append('g')
    .attr('transform', `translate(${2 * padding},${textHeight + padding})`);

  barPlot.append('text')
    .attr('fill', 'lightgrey')
    .attr('font-size', `${textLabelHeight}px`)
    .text('Ø L2-distance to µ')
    .attr('opacity', 0)
    .transition()
      .delay(500)
      .attr('opacity', 1)

  const x = scaleBand()
    .domain(groups)
    .range([0, barWidth])
    .padding(0.2);

  const y = scaleLinear()
    .domain([0, maxDist])
    .range([barHeight, 0]);

  const xAxis = axisBottom(x)
    .tickSizeInner(innerTickSize)
    .tickSizeOuter(outerTickSize)
    .tickPadding(tickPadding);

  barPlot.append('g')
    .attr('transform', `translate(0,${textLabelHeight + barHeight})`)
    .attr('class', 'axis')
    .attr('opacity', 0)
    .transition()
    .delay(500)
      .attr('opacity', 1)
    .call(xAxis);

  const yAxis = axisLeft(y)
    .tickFormat(d => labelFormatter(d))
    .tickSizeInner(innerTickSize)
    .ticks(barHeight / tickSize)
    .tickSizeOuter(outerTickSize)
    .tickPadding(tickPadding);

  barPlot.append('g')
    .attr('transform', `translate(0,${textLabelHeight})`)
    .attr('class', 'axis')
    .attr('opacity', 0)
    .transition()
      .delay(500)
      .attr('opacity', 1)
    .call(yAxis);

  const bars = barPlot.append('g')

  bars.selectAll('rect')
    .data(barData)
    .join('rect')
      .attr('x', (d, i) => x(groups[i]))
      .attr('y', textLabelHeight + barHeight)
      .attr('width', x.bandwidth())
      .attr('height', 0)
      .attr('fill', (d, i) => i === 0 ? cluster.color : 'grey')
      .transition()
      .delay(500)
      .duration(1000)
        .attr('y', d => textLabelHeight + y(d))
        .attr('height', d => barHeight - y(d))
};