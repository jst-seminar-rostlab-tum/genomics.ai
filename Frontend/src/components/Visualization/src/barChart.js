import * as d3 from "d3";
import "./bar.css"

export const addBarPlot = (BarContainer, data) => {
  const svg = d3.select(BarContainer).append("svg");
  const margin = 80;
  const w = 300;
  const h = 300;

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
  const info = Array.from(groupedByCelltypes).filter(d => d[0] !== undefined).map(d => [d[0], d[1].length]);
  const labels = Array.from(groupedByCelltypes.keys()).filter(d => d !== undefined);
  const xScale = d3.scaleBand().domain(labels).range([margin, w]).padding(0.4);
  const yScale = d3.scaleLinear().domain([0, data.length]).range([h, margin + 25]);


  const g = svg.append("g")

  //Titel of grafic
  g.append("text")
    .attr("class", "title")
    .attr("x", margin + 40)
    .attr("y", 15)
    .attr("font-size", "15px")
    .text("#Cells per Category")
    .style("font-family", "inherit")
    .style("font-weight", 640)

  //Axes
  const xAxis = d3.axisBottom(xScale)
  // .ticks([labels]);
  const yAxis = d3.axisLeft(yScale).tickFormat(d => d).ticks(5).tickSize(-220);


  g.append("g")
    .attr("id", "x-axis")
    .attr("class", "x-axis")
    .attr("transform", "translate(0, " + (h - margin) + ")")
    .call(xAxis)
    .selectAll("text")
    .attr("transform", "translate(-10,0)rotate(-45)")
    .style("text-anchor", "end")
    .style("color", "black")


  g.append("g")
    .attr("id", "y-axis")
    .attr("class", "y-axis")
    .attr("transform", "translate(" + margin + "," + (-margin) + ")")
    .call(yAxis)
    .selectAll("text")
    .style("text-anchor", "end")
    .style("color", "black")


  g.append("g")
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
      return h - yScale(d[1])
    });

  d3.select("#bars")
    .attr("transform", "translate( 0, " + (-margin) + ")")
}