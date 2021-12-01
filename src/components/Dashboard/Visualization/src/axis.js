import * as d3 from "d3";
import { width, height } from "./constants";

const tickSize = 50;
const tickPadding = 5;

const customizeAxis = (axis, length, gridLineLength) => {
  axis
    .tickSizeInner(-gridLineLength)
    .ticks(length / tickSize)
    .tickPadding(tickPadding)
    .tickSizeOuter(0);
}

export const addAxisLabel = (selection, text, x, y, rotateDeg = 0) => {
  selection
    .append("text")
    .attr("transform", `rotate(-${rotateDeg})`)
    .attr("class", "axislabel")
    .attr("x", x)
    .attr("y", y)
    .text(text)
}

export const addAxes = (gX, gY, domX, domY) => {
  const xAxis = d3.axisBottom(domX);
  const yAxis = d3.axisRight(domY);
  customizeAxis(xAxis, width, height);
  customizeAxis(yAxis, height, width);
  gX.call(xAxis);
  gY.call(yAxis);
  return { xAxis, yAxis };
}