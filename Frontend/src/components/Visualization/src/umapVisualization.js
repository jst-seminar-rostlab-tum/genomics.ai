import * as d3 from "d3";
import * as cons from "./constants";
import { zoomM } from "./newZoom";
import "./tooltip.css";
import { addBarPlotCell, addBarPlotBatch, addBarPlotPredict } from "./barChart"
import {getColoringModes, setColoring} from "./coloring"

//Create a color scale based on the chosen mode and its values

const getMin = (w, h) => {
  return d3.min([h, w]);
}

const addGroup = (svg, id) => {
  return svg.append('g').attr('id', `${id}`);
};

export class UmapVisualization2 {

  constructor(container, data) {
    d3.select(container).selectAll("*").remove();
    this.svg = d3.select(container).append('svg')
        .attr('id', 'vis_svg');
    this.gCells = addGroup(this.svg, 'cells');
    // this.gLabels = addGroup(this.svg, 'labels');
    this.coloringModes = getColoringModes(data);
    this.tooltip = d3.select(container).append("div");
    this.mode = undefined;
    this.data = data;
    this.graphs = ["batch", "cell"];
    this.hiddenCells = [];
  };

  getAvailableGraphs() {
    return this.graphs;
  }

  drawGraph(container, graph, width, height) {
    switch (graph) {
      case "batch":
        if(Object.keys(this.coloringModes).includes("batch")){
          d3.select(container).append("div");
          addBarPlotBatch(container, this.data, width, height);
        }
        break;
      case "cell":
        if(Object.keys(this.coloringModes).includes("cell_type")){
          d3.select(container).append("div");
          this.barChartCell = addBarPlotCell(container, this.data, width, height);
        }
        break;
      default:
        break;
    }
  };

  isHidden(cell){
    for (let i = 0; i < this.hiddenCells.length; i++){
       if (cell[this.hiddenCells[i][0]] === this.hiddenCells[i][1]) {
         return true;
       }
    }
    return false;
  }

  //Hide
  after(category, value) {
    this.addHiddenCell(category, value);
    this.hideShowCells();
  }

  //Show
  before(category, value) {
    this.deleteHiddenCell(category, value);
    console.log(this.hiddenCells);
    this.hideShowCells();
  }

  beforeAll() {
    this.hiddenCells = [];
    this.cells
      .style("visibility", "visible");
  }

  hideShowCells() {
    this.cells
      .style("visibility", (d) => { return this.isHidden(d) ? "hidden" : "visible" });
  }

  addHiddenCell(category, value) {
    this.hiddenCells.push([category, value]);
  }

  filterCells(cell) {
    cell[0] == this[0] && cell[1] == this[1];
  }

  deleteHiddenCell(category, value) {
    this.hiddenCells = this.hiddenCells.filter(function (cell) {
      return !(cell[0] == category && cell[1] == value)});
  }

  //////////////////////////////////////////////////// These should be deleted and the use in F3 changed, when
                                                      // it's certain how the attributes are going to be named
  showReference() {
    this.before("type", "reference");
  }

  showQuery() {
    this.before("type", "query");
  }

  hideQuery() {
    this.after("type", "query");
  }

  hideReference() {
    this.after("type", "reference");
  }
////////////////////////////////////////////////////////
  predictedCellsTransparent() {
    this.cells
      .style("opacity", (d) => {return d["predictions"] === "True" ? 0.3 : 0.85 });
  }

  predictedCellsVisible() {
    this.cells
      .style("opacity", 1);
  }

  //Set a color mode
  setColorMode(mode) {
    const colorScale = setColoring(mode, this.data);
    if (this.cells != null) {

      this.cells.style("fill", (d) => colorScale(d[mode]));
    }
    this.mode = mode;
    this.colorScale = colorScale;
  }

  //Sizing and resizing the cells
  resize(w, h) {
    if(this.cells === undefined) {
      return this.render(w, h);
    }

    const min = getMin(w, h);
    const data = this.data;
    const r = 0.003 * min;

    const xScaleHelper = d3.scaleLinear()
      .domain([d3.min(data.map(d => parseFloat(d.x))), d3.max(data.map(d => parseFloat(d.x)))])
      .range([cons.margin, min - cons.margin]);

    const translate = (w - cons.margin - (xScaleHelper(d3.max(data.map(d => parseFloat(d.x)))))) / 2;

    //Scales
    const xScale = d3.scaleLinear()
      .domain([d3.min(data.map(d => parseFloat(d.x))), d3.max(data.map(d => parseFloat(d.x)))])
      .range([translate, min - cons.margin + translate]);

    const yScale = d3.scaleLinear()
      .domain([d3.min(data.map(d => parseFloat(d.y))), d3.max(data.map(d => parseFloat(d.y)))])
      .range([cons.margin, min - cons.margin]);

    this.cells
      .attr("cx", d => xScale(parseFloat(d.x)))
      .attr("cy", d => yScale(parseFloat(d.y)))

      this.cells
      .attr("r", r);


    return [xScale, yScale];
  }

  //Constructing the svg
  async render(w, h) {
    const data = this.data;

    const min = getMin(w, h);

    const r = 0.003 * min;

    const _this = this;

    //svg
    this.svg
      .attr("width", w)
      .attr("height", h)

    //Circle cell
    this.gCells.selectAll("*").remove();

    this.cells =
      this.gCells
        .selectAll("cell")
        .data(data)
        .enter()
        .append("circle")
        .attr("r", r)
        .style("fill", (d) => (_this.mode ? _this.colorScale(d[_this.mode]) : "black"))
        .style("opacity", cons.originalOpacity)


    //Setting coordinates
    const [xScale, yScale] = this.resize(w, h);

    //tooltip
    var tooltip =
      this.tooltip
        .attr('class', 'tooltip')
        .attr('opacity', 1)
        .style("visibility", "hidden")


    this.cells
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
      })


    //Pan and mouse zoom
    this.svg
      .attr('class', 'mouse-capture')
      .attr('x', -w)
      .attr('y', -h)
      .attr('width', w)
      .attr('height', h)
      .style('fill', 'white')
      .lower()
      .call(zoomM)

    this.cells
      .attr("cx", d => xScale(parseFloat(d.x)))
      .attr("cy", d => yScale(parseFloat(d.y)))

  }

}
