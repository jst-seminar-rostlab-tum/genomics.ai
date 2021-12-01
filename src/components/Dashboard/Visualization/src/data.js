import { width, height } from "./constants";
import * as d3 from "d3";

const displayedPoints = 2500;

export const prepareData = ([metaData, umapData]) => {
  const clusterMapping = new Map();

  metaData
    .filter(({ CELLTYPE }) => isNaN(parseFloat(CELLTYPE)))
    .forEach(m => clusterMapping.set(parseInt(m.CLUSTER), m.CELLTYPE));

  const pointsData = umapData
    .slice(0, displayedPoints)
    .map(d => ({
    id: d.cellType,
    x: d.x,
    y: d.y,
    cluster: parseInt(d.cluster)
    }))
    .filter(d => clusterMapping.has(d.cluster));

  return { clusterMapping, pointsData };
}

const fraction = 0.1;

export const getDomains = (clusterMapping, pointsData) => {
  const clusterLabels = clusterMapping.keys();

  const xMin = d3.min(pointsData, d => +d.x);
  const xMax = d3.max(pointsData, d => +d.x);
  const yMin = d3.min(pointsData, d => +d.y);
  const yMax = d3.max(pointsData, d => +d.y);

  const lagX = Math.ceil(fraction * (xMax - xMin));
  const lagY = Math.ceil(fraction * (yMax - yMin));

  const domX = d3
    .scaleLinear()
    .domain([xMin - lagX, xMax + lagX])
    .range([0, width]);

  const domY = d3
    .scaleLinear()
    .domain([yMin - lagY, yMax + lagY])
    .range([height, 0]);

  const colorScale = d3
    .scaleOrdinal()
    .domain(clusterLabels)
    .range(d3.schemeCategory10);

  return { domX, domY, colorScale };
}