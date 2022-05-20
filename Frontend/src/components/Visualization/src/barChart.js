import * as d3 from "d3";
import "./bar.css"

const groupBy = (data, cat) => {
  let list = []
  let groupedBy = new Map();
  for (const row in data) {
    const item = data[row];

    if (groupedBy.has(item[cat])) {
      list = groupedBy.get(item[cat]);
      list.push(item);
      groupedBy.set(item[cat], list);

    } else {
      groupedBy.set(item[cat], [item]);
    }
  }

  return groupedBy;
}


const addBarPlot = (barContainer, data, groupedBy, title, w, h) => {
  const svg = d3.select(barContainer).append("svg");
  const marginBottom = 85;
  const plotTitleOffset = 25;
  const plotTitleCentered = 25;
  //const w = 270;
  //const h = 270;

  svg
    .attr("width", w)
    .attr("height", h)

  //Scales
  const info = Array.from(groupedBy).filter(d => d[0] !== undefined).map(d => [d[0], d[1].length]);
  const labels = Array.from(groupedBy.keys()).filter(d => d !== undefined);
  const xScale = d3.scaleBand().domain(labels).range([marginBottom, w]).padding(0.4);
  const yScale = d3.scaleLinear().domain([0, data.length]).range([h, marginBottom + plotTitleOffset]);


  const g = svg.append("g")
  let fontSizeTitle;
  let fontSizeOther;
  if (h > 350) {
    fontSizeTitle = "20px";
    fontSizeOther = "15px";
  } else {
    fontSizeTitle = "15px";
    fontSizeOther = "10px";
  }
  //Titel of grafic
  g.append("text")
    .attr("class", "title")
    .attr("x", marginBottom + plotTitleCentered)
    .attr("y", 15)
    .attr("font-size", fontSizeTitle)
    .text(title)
    .style("font-family", "inherit")
    .style("font-weight", 640)

  //Axes
  const xAxis = d3.axisBottom(xScale);
  const yAxis = d3.axisLeft(yScale).tickFormat(d => d).ticks(5).tickSize(-h);


  g.append("g")
    .attr("id", "x-axis")
    .attr("class", "x-axis")
    .attr("transform", "translate(0, " + (h - marginBottom) + ")")
    .call(xAxis)
    .selectAll("text")
    .attr("transform", "translate(-10,0)rotate(-45)")
    .style("text-anchor", "end")
    .style("color", "black")
    .attr("font-size", fontSizeOther);


  g.append("g")
    .attr("id", "y-axis")
    .attr("class", "y-axis")
    .attr("transform", "translate(" + marginBottom+ "," + (-marginBottom) + ")")
    .call(yAxis)
    .selectAll("text")
    .style("text-anchor", "end")
    .style("color", "black")
    .attr("font-size", fontSizeOther);


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
    })
    .attr("transform", "translate( 0, " + (-marginBottom) + ")")  

    g.append("g")
    .attr("id", "labels")
    .selectAll("label")        
    .data(info)
    .enter()
    .append("text")
    .attr("class","label")
    .attr("x", (d) => xScale(d[0]))
    .attr("y", (d) => yScale(d[1]) - marginBottom - 13)
    .attr("dy", ".75em")
    .text((d) => d[1])
    .attr("fill", "black")
    .attr("font-size", fontSizeOther);
     
}

export const addBarPlotCell = (barContainer, data, w, h) =>{
  const groupedByCellTypes = groupBy(data, "cell_type");
  addBarPlot(barContainer, data, groupedByCellTypes, "#Cells per Cell Type", w, h);
}

export const addBarPlotBatch = (barContainer, data, w, h) =>{
  const groupedByBatch = groupBy(data, "batch");
  addBarPlot(barContainer, data, groupedByBatch, "#Cells per Batch", w, h);
}