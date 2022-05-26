import * as d3 from "d3";
import "./bar.css"
import * as cons from "./constants";


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

  // Defines the scales
  const info = Array.from(groupedBy).filter(d => d[0] !== undefined).map(d => [d[0], d[1].length]);
  const labels = Array.from(groupedBy.keys()).filter(d => d !== undefined);
  const xScale = d3.scaleBand().domain(labels).range([marginBottom, w]).padding(0.4);
  const yScale = d3.scaleLinear().domain([0, data.length]).range([h, marginBottom + cons.plotTitleOffset]);


  const g = svg.append("g")

  // Titel of grafic
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

  /*bars
      .on("mouseover", function (m, d) { //getter bigger on hover
        d3.select(this)
          .transition()
          .duration(100)
          .attr('r', r * 1.5)
          .style("stroke", "black")
          .style("opacity", 1)
        if (!_this.mode) return;
        const att = d[_this.mode];
        const xPos = parseFloat(m.x) + 5;
        const yPos = parseFloat(m.y) - 35;

        tooltip
          .style("visibility", "visible")

        tooltip
          .html(att)
          .style("top", `${yPos}px`)
          .style("left", `${xPos}px`)

      })
      .on("mouseout", function () {
        tooltip.style("visibility", "hidden");
        d3.select(this)
          .transition()
          .duration(100)
          .attr('r', r)
          .style("stroke", "none")
          .style("opacity", 0.8)
      })*/

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