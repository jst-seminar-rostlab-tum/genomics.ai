// const { margin } = require("@mui/system");
import * as d3 from "d3";
import "./bar.css"

export const addBarPlot = (BarContainer, data) => {
  const svg = d3.select(BarContainer).append("svg");
  const margin = 10;
  const w = 250;
  const h = 250;

  svg
    .attr("width", w)
    .attr("height", h)

  let groupedByCelltypes = new Map();
  let groupedByBatch = new Map();

  let list = []
  for (const row in data) {
    const item = data[row];

    if (groupedByCelltypes.has(item.cell_type)) {
      list = groupedByCelltypes.get(item.cell_type);
      list.push(item);
      groupedByCelltypes.set(item.cell_type, list);

    } else {
      groupedByCelltypes.set(item.cell_type, [item]);
    }
  }

  for (const row in data) {
    const item = data[row];

    if (groupedByBatch.has(item.batch)) {
      list = groupedByBatch.get(item.batch);
      list.push(item);
      groupedByBatch.set(item.batch, list);

    } else {
      groupedByBatch.set(item.batch, [item]);
    }
  }

  //Scales
  const info = Array.from(groupedByCelltypes).filter(d => d[0] !== undefined).map( d => [d[0],d[1].length]);
  const labels = Array.from(groupedByCelltypes.keys()).filter( d => d!==undefined);
  const xScale = d3.scaleBand().domain(labels).range([margin, w - margin]);
  const yScale = d3.scaleLinear().domain([0, data.length]).range([h - margin, 0]);

  //Axes
  const xAxis = d3.axisBottom(xScale).ticks([labels]);
  const yAxis = d3.axisLeft(yScale);

  svg.append("g")
    .attr("id", "y-axis")
    .attr("class", "tick")
    .attr("transform", "translate(" + margin + ", 0)")
    .call(yAxis);

  svg.append("g")
    .attr("id", "x-axis")
    .attr("class", "tick")
    .attr("transform", "translate( 0, " + (h - margin) + ")")
    .call(xAxis);

  svg
  .append("g")
  .attr("id", "bars")
  .selectAll("rect")
    .data(info)
    .enter()
    .append("rect")
    .attr("class", "bar")
    .attr("x", (d) => xScale(d[0]))
    .attr("y", (d) => yScale(d[1]))
    .attr("width", xScale.bandwidth())
    .attr("height", (d) => {
      console.log(h - yScale(d[1]));
      return h - yScale(d[1])})
    .append('title')
    .text(d => {console.log(d[1]/data.length);
      return d[1]/data.length;})

    d3.select("#bars")
    .attr("transform", "translate( 0, " + (-margin)  + ")");

  // const xScale = d3.scaleBand().range([margin, w - margin]).padding(0.4);
  // xScale.domain(data.map(function (d) {
  //   return d.cell_type;
  // }));

  // const yScale = d3.scaleLinear().range([h - margin,0]);
  // yScale.domain([0, data.length]);

  // g.append("g")
  //   .attr("transform", "translate(0," + (h-2*margin) + ")")
  //   .call(d3.axisBottom(xScale))
  //   .selectAll("text")
  //   .attr("transform", "translate(-10,0)rotate(-45)")
  //   .style("text-anchor", "end")


  // g.append("g")
  //   .call(d3.axisLeft(yScale).tickFormat(function (d) {
  //     return d;
  //   })
  //     .ticks(3))
  //   .append("text")
  //   .attr("transform", "rotate(-90)")
  //   .attr("y", 6)
  //   .attr("dy", "-5.1em")
  // // .attr("text-anchor", "end")
  // // .attr("stroke", "black")
  // // .text("#All Cells");

  // console.log(groupedByCelltypes);

  // g.selectAll("rect")
  //   .data(data)
  //   .enter()
  //   .append("rect")
  //   .attr("class", "bar")
  //   .attr("x", function (d) {
  //     return xScale(d.cell_type);
  //   })
  //   .attr("y", function (d) {

  //     //TODO: switch between cells Type and Batch
  //     return yScale(getNum(d.cell_type));
  //   })
  //   .attr("width", xScale.bandwidth())
  //   .attr("height", function (d) {
  //     // return height - yScale(getNum(d.cell_type));
  //     return yScale(getNum(d.cell_type));
  //   })
  //   // .attr("fill", "steelblue")
  //   .append('title')
  //   .text(d => getNum(d.cell_type) / data.length)


}
//   );
// }
