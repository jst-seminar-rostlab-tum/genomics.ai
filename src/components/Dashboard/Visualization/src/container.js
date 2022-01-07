import { margin, width, height } from './constants';

const addBorderLine = (selection, x1, y1, x2, y2) => {
  selection
    .append('line')
    .style('stroke', 'black')
    .attr('x1', x1)
    .attr('y1', y1)
    .attr('x2', x2)
    .attr('y2', y2);
};

const addClipPath = (selection) => {
  selection
    .append('defs')
    .append('svg:clipPath')
    .attr('id', 'clip')
    .append('svg:rect')
    .attr('width', width)
    .attr('height', height);
};

const legendWidth = 200;

export const prepareContainer = (container) => {
  const svgDots = container
    .append('svg')
    .attr('width', margin.left + width + margin.right)
    .attr('height', margin.top + height + margin.bottom);

  addClipPath(svgDots);
  addBorderLine(svgDots, margin.left, margin.top, margin.left + width, margin.top);
  addBorderLine(svgDots, margin.left, margin.top, margin.left, margin.top + height);

  const svgLegend = container
    .append('svg')
    .attr('width', legendWidth)
    .attr('height', margin.top + height + margin.bottom);

  const buttons = container
    .append('div')
    .attr('id', 'buttons');

  const tooltip = container
    .append('div')
    .attr('class', 'tooltip')
    .style('display', 'none');

  return {
    svgDots, svgLegend, buttons, tooltip,
  };
};

export const addGroups = (svgDots) => {
  const gX = svgDots.append('g')
    .attr('class', 'axis axis--x')
    .attr('transform', `translate(${margin.left},${height + margin.top})`);

  const gY = svgDots.append('g')
    .attr('class', 'axis axis--y')
    .attr('transform', `translate(${margin.left + width},${margin.top})`);

  const gDots = svgDots.append('g')
    .attr('clip-path', 'url(#clip)')
    .attr('transform', `translate(${margin.left},${margin.top})`);

  return { gX, gY, gDots };
};
