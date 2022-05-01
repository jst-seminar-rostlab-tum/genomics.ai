// const { margin } = require("@mui/system");
import * as d3 from "d3";

export const addBarPlot = (BarContainer, data) => {
  const svg = d3.select(BarContainer).append("svg");
    const margin = 0;
    const w = 250;
    const h = 250;

    svg
    .attr("width", w)
    .attr("height", h)

  // svg.append("text")
  //   .attr("transform", "translate(100,0)")
  //   .attr("x", 0)
  //   .attr("y", 50)
  //   .attr("font-size", "20px")
  //   .text("#Gene per Percentile")

  const xScale = d3.scaleBand().range([0, w]) //.padding(0.4),
    const yScale = d3.scaleLinear().range([0, h]);

  // const g = svg.append("g")
    // .attr("transform", "translate(" + 100 + "," + 100 + ")");

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

    const getNum = (type) => {
      // if (mode === cell_type)
      return [...groupedByCelltypes.get(type)].length;
      // if (Mode === batch)
      //   return [...groupedByBatch.get(type)].length;
    }

    const getBatchNum = (type) => {
      // if (mode === cell_type)
      // return [...groupedByCelltypes.get(type)].length;
      // if (Mode === batch)
      return [...groupedByBatch.get(type)].length;
    }

    xScale.domain(data.map(function (d) {
      return d.cell_type;
    }));
    yScale.domain([0, data.length]);

    g.append("g")
      // .attr("transform", "translate(0," + height + ")")
      .call(d3.axisBottom(xScale))
      .selectAll("text")
      .attr("transform", "translate(-10,0)rotate(-45)")
      .style("text-anchor", "end")


    g.append("g")
      .call(d3.axisLeft(yScale).tickFormat(function (d) {
        return d;
      })
        .ticks(3))
      .append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", 6)
      .attr("dy", "-5.1em")
    // .attr("text-anchor", "end")
    // .attr("stroke", "black")
    // .text("#All Cells");

    g.selectAll("rect")
      .data(data)
      .enter()
      .append("rect")
      .attr("class", "bar")
      .attr("x", function (d) {
        return xScale(d.cell_type);
      })
      .attr("y", function (d) {

        //TODO: switch between cells Type and Batch
        return yScale(getNum(d.cell_type));
      })
      .attr("width", xScale.bandwidth())
      .attr("height", function (d) {
        // return height - yScale(getNum(d.cell_type));
        return yScale(getNum(d.cell_type));
      })
      .attr("fill", "steelblue")
      .append('title')
      .text(d => getNum(d.cell_type) / data.length)


  }
//   );
// }
