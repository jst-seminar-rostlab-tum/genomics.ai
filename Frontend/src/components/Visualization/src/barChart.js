import * as d3 from "d3";
import "./bar.css"
import * as cons from "./constants";

/*
Most of the styling is in the bar.css file
*/

// Groups the data based on a category specified in cat
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

//Randomizes the row of the bars/values, to reduce label overlapping
const shuffleArray = (array) => {
  for (var i = array.length - 1; i > 0; i--) {
      var j = Math.floor(Math.random() * (i + 1));
      var temp = array[i];
      array[i] = array[j];
      array[j] = temp;
  }
}

// Constructs the graph
const addBarPlot = (barContainer, data, groupedBy, title, cWidth, cHeight) => {
  const svg = d3.select(barContainer).append("svg");

  let h = cHeight;
  let w = cWidth;
  svg
    .attr("width", w)
    .attr("height", h)

  // Changes the font size/bottom margin based on the size of the graph (mostly important for the pop up)
  let fontSizeTitle;
  let fontSizeOther;
  let marginBottom;
  let big = h > 350; // to change the pop up based on if it's big or not
  if (big) {
    fontSizeTitle = cons.fontSizeTBig;
    fontSizeOther = cons.fontSizeOBig;
    marginBottom = cons.marginBottomBig;
  } else {
    fontSizeTitle = cons.fontSizeTSmall;
    fontSizeOther = cons.fontSizeOSmall;
    marginBottom = cons.marginBottomSmall;
  }
  
  const numberOfGraphValues = 15; // Has to be one more than specified

  const sortedData = Array.from(groupedBy).sort((a, b) => {if (a[1] === b[1]) {
        return 0;
    }
    else {
        return (a[1] > b[1]) ? -1 : 1;
    }});
  const splicedData = sortedData.slice(0, sortedData.length >= numberOfGraphValues ? numberOfGraphValues : sortedData.length);
  shuffleArray(splicedData);

  // Defines the scales
  const info = splicedData.filter(d => d[0] !== undefined).map(d => [d[0], d[1].length]);

  const max = d3.max(info.map(d => d[1])) * 1.1;

  const labels = splicedData.map(([k, v]) => k).filter(d => d !== undefined); //gets the names of the category values

  //Scales
  const xScale = d3.scaleBand().domain(labels).range([marginBottom, w]).padding(0.4);
  const yScale = d3.scaleLinear().domain([0, max]).range([h, marginBottom + cons.plotTitleOffset]);

  const g = svg.append("g")

  //Title of the graph
  g.append("text")
    .attr("class", "title")
    .attr("x", marginBottom + cons.plotTitleCentered)
    .attr("y", 15)
    .attr("font-size", fontSizeTitle)
    .text(title)
    .style("font-family", "inherit")
    .style("font-weight", 640)

  // Axes
  const xAxis = d3.axisBottom(xScale);
  const yAxis = d3.axisLeft(yScale).tickFormat(d => d).ticks(5).tickSize(-h);

  // Makes the x-axis
  g.append("g")
    .attr("id", "x-axis")
    .attr("class", "x-axis")
    .attr("transform", "translate(0, " + (h - marginBottom) + ")")
    .call(xAxis)
    .selectAll("text")
    //styling ticks
    .attr("transform", "translate(-10,0)rotate(-45)")
    .style("text-anchor", "end")
    .style("color", "black")
    .attr("font-size", fontSizeOther)

  // Makes the y-axis
  g.append("g")
    .attr("id", "y-axis")
    .attr("class", "y-axis")
    .attr("transform", "translate(" + marginBottom + "," + (-marginBottom) + ")")
    .call(yAxis)
    //styling ticks
    .selectAll("text")
    .style("text-anchor", "end")
    .style("color", "black")
    .attr("font-size", fontSizeOther);

  // Draws the bars
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
    .append("title")
    .text((d) => (d[1]))

  // Only the big pop up has labels
  if (big) {
    // Adds the values on top of the bars
    g.append("g")
      .attr("id", "labels")
      .selectAll("label")
      .data(info)
      .enter()
      .append("text")
      .attr("class", "label")
      .attr("x", (d) => xScale(d[0]))
      .attr("y", (d) => yScale(d[1]) - marginBottom - 13)
      .attr("dy", ".75em")
      .text((d) => d[1])
      .attr("fill", "black")
      .attr("font-size", fontSizeOther)
  }

}

// Draws the cell_type chart
export const addBarPlotCell = (barContainer, data, w, h) => {
  const groupedByCellTypes = groupBy(data, "cell_type");
  addBarPlot(barContainer, data, groupedByCellTypes, "#Cells per Cell Type", w, h);
}

// Draws the batch chart
export const addBarPlotBatch = (barContainer, data, w, h) => {
  const groupedByBatch = groupBy(data, "batch");
  addBarPlot(barContainer, data, groupedByBatch, "#Cells per Batch", w, h);
}