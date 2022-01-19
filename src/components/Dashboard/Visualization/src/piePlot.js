import { pie, interpolate, arc } from 'd3';
import {
  innerPieRadius, outerPieRadius, padding, barWidth, textHeight, textLabelHeight,
} from './constants';

const sampleLabel = '# Samples';

export const addPiePlot = (labelContainer, cluster, clusterCollection) => {
  const generator = pie().value((d) => d.count);
  const pieData = generator(clusterCollection.clusters);
  const angleInterpolation = interpolate(generator.startAngle()(), generator.endAngle()());
  const innerRadiusInterpolation = interpolate(0, innerPieRadius);
  const outerRadiusInterpolation = interpolate(0, outerPieRadius);

  const piePlot = labelContainer.append('g')
    .attr('transform', `translate(${3 * padding + barWidth},${textHeight + padding})`);

  piePlot.append('text')
    .attr('x', outerPieRadius - 7 * sampleLabel.length / 2)
    .attr('fill', 'lightgrey')
    .attr('font-size', `${textLabelHeight}px`)
    .text(sampleLabel)
    .attr('opacity', 0)
    .transition()
    .delay(500)
    .attr('opacity', 1);

  const countString = `${cluster.count}`;

  piePlot.append('text')
    .attr('x', outerPieRadius - 7 * countString.length / 2)
    .attr('y', 54)
    .attr('fill', 'lightgrey')
    .attr('font-size', `${textLabelHeight}px`)
    .text(`${cluster.count}`)
    .attr('opacity', 0)
    .transition()
    .delay(1000)
    .attr('opacity', 1);

  const { totalCount } = clusterCollection;
  const percentage = 100 * cluster.count / totalCount;
  const percentageString = `${percentage.toFixed(1)}%`;

  piePlot.append('text')
    .attr('x', outerPieRadius - 7 * percentageString.length / 2)
    .attr('y', 71)
    .attr('fill', 'lightgrey')
    .attr('font-size', `${textLabelHeight}px`)
    .text(percentageString)
    .attr('opacity', 0)
    .transition()
    .delay(1000)
    .attr('opacity', 1);

  const chart = piePlot.append('g')
    .attr('transform', `translate(${outerPieRadius},${outerPieRadius + textLabelHeight})`);

  const arcs = chart.selectAll('path')
    .data(pieData)
    .enter()
    .append('path')
    .attr('fill', (d, i) => (clusterCollection.clusters[i].name === cluster.name ? cluster.color : 'lightgrey'));

  const ar = arc();

  arcs.transition()
    .duration(1500)
    .attrTween('d', (d) => {
      const originalEnd = d.endAngle;
      return (t) => {
        const currentAngle = angleInterpolation(t);
        if (currentAngle < d.startAngle) return '';
        d.endAngle = Math.min(currentAngle, originalEnd);
        return ar(d);
      };
    });

  chart
    .transition()
    .duration(1500)
    .tween('arcRadii', () => (t) => ar
      .innerRadius(innerRadiusInterpolation(t))
      .outerRadius(outerRadiusInterpolation(t)));
};
