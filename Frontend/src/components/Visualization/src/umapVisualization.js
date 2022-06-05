import * as d3 from "d3";
import * as cons from "./constants";
import { zoomM } from "./newZoom";
import "./tooltip.css";
import { addBarPlotCell, addBarPlotBatch, addBarPlotPredict } from "./barChart"
import { getColoringModes, setColoring } from "./coloring"

//Gets the minimum out of 2 values - used to determine the smaller side of the container, to make scatter plot square
const getMin = (w, h) => {
  return d3.min([h, w]);
}

//Adds a group
const addGroup = (svg, id) => {
  return svg.append('g').attr('id', `${id}`);
};

export class UmapVisualization2 {

  constructor(container, data) {
    d3.select(container).selectAll("*").remove();
    this.svg = d3.select(container).append('svg')
      .attr('id', 'vis_svg');
    this.gCells = addGroup(this.svg, 'cells');
    this.coloringModes = getColoringModes(data);
    this.tooltip = d3.select(container).append("div");
    this.mode = undefined; //coloring mode
    this.data = data;
    this.graphs = ["batch", "cell"]; //constant attributes expected cells to have
    this.hiddenCells = [];
  };

  // Returns all available graphs
  getAvailableGraphs() {
    //If the cells don't have a cell type, as with totalVI, we only have a batch graph
    return (Object.keys(this.coloringModes).includes("cell_type"))? this.graphs: ["batch"];
  }

  // Draws a graph, which is in the this.graphs list 
  drawGraph(container, graph, width, height) {
    switch (graph) {
      case this.graphs[0]: // Draws the batch graph
        if (Object.keys(this.coloringModes).includes("batch")) {
          d3.select(container).append("div");
          addBarPlotBatch(container, this.data, width, height);
        }
        break;
      case this.graphs[1]: // Draws the cell type graph
        if (Object.keys(this.coloringModes).includes("cell_type")) {
          d3.select(container).append("div");
          this.barChartCell = addBarPlotCell(container, this.data, width, height);
        }
        break;
      default:
        break;
    }
  };

  // Checks if a cell fullfills one of the [category, value] of the this.hiddenCells list
  // If yes -> true
  // If not -> false
  isHidden(cell) {
    for (let i = 0; i < this.hiddenCells.length; i++) {
      if (cell[this.hiddenCells[i][0]] === this.hiddenCells[i][1]) {
        return true;
      }
    }
    return false;
  }

  // Hide the cells that have the given value in the given category
  after(category, value) {
    this.addHiddenCell(category, value);
    this.hideShowCells();
  }

  // Show the cells that have the given value in the given category
  before(category, value) {
    this.deleteHiddenCell(category, value);
    this.hideShowCells();
  }

  // Resets the visibility for all cells
  beforeAll() {
    this.hiddenCells = [];
    this.cells
      .style("visibility", "visible");
  }

  // Changes the visibility of the cells based on the content of this.hiddenCells
  hideShowCells() {
    this.cells
      .style("visibility", (d) => { return this.isHidden(d) ? "hidden" : "visible" });
  }

  // Add a [category, value] to the hidden cells
  addHiddenCell(category, value) {
    this.hiddenCells.push([category, value]);
  }

  filterCells(cell) {
    cell[0] == this[0] && cell[1] == this[1];
  }

  // Removes a [category, value] from the hidden cells
  deleteHiddenCell(category, value) {
    this.hiddenCells = this.hiddenCells.filter(function (cell) {
      return !(cell[0] == category && cell[1] == value)
    });
  }

  // Unhides the reference cells
  showReference() {
    this.before("type", "reference");
  }

  // Unhides the query cells
  showQuery() {
    this.before("type", "query");
  }

  // Hides the query cells
  hideQuery() {
    this.after("type", "query");
  }

  // Hides the reference cells
  hideReference() {
    this.after("type", "reference");
  }

  /* The following methods are their own methods to allow for 
  the option to add the possibility of changing the opacity based 
  on multiple attributes */

  // Sets the opacity of all predicted cells to a reduced value
  predictedCellsTransparent() {
    this.reduceOpacity("predicted", "yes");
  }

  // Sets the opacity of all predicted cells to the original value
  predictedCellsVisible() {
    this.resetOpacity();
  }

  // Reduces the opacity of the cells that have the given value in the given category
  reduceOpacity(category, value) {
    this.cells
      .style("opacity", (d) => { return d[category] === value ? cons.reducedOpacity : cons.originalOpacity });
  }

  // Resets the opacity of all cells (back to the original opacity)
  resetOpacity() {
    this.cells
      .style("opacity", cons.originalOpacity);
  }

  // Colors the cells based on the specified mode
  setColorMode(mode) {
    const colorScale = setColoring(mode, this.data);
    if (this.cells != null) {

      this.cells.style("fill", (d) => colorScale(d[mode]));
    }
    this.mode = mode;
    this.colorScale = colorScale;
  }

  // Sizing and resizing the cells
  resize(w, h) {
    if (this.cells === undefined) {
      return this.render(w, h);
    }

    const min = getMin(w, h);
    const data = this.data;
    const r = 0.003 * min;

    //Makes an x-axis scale helper, using min and max value as domain, and the distance from margin to margin for range 
    const xScaleHelper = d3.scaleLinear()
      .domain([d3.min(data.map(d => parseFloat(d.x))), d3.max(data.map(d => parseFloat(d.x)))])
      .range([cons.margin, min - cons.margin]);

    //As the previous scale makes the scatterplot off-center, we compute a translation distance
    const translate = (w - cons.margin - (xScaleHelper(d3.max(data.map(d => parseFloat(d.x)))))) / 2;

    // Scales
    //Official x-axis scale, adding the transaltion distance
    const xScale = d3.scaleLinear()
      .domain([d3.min(data.map(d => parseFloat(d.x))), d3.max(data.map(d => parseFloat(d.x)))])
      .range([translate, min - cons.margin + translate]);

    //y-axis scale
    const yScale = d3.scaleLinear()
      .domain([d3.min(data.map(d => parseFloat(d.y))), d3.max(data.map(d => parseFloat(d.y)))])
      .range([cons.margin, min - cons.margin]);

    //Setting attributes
    this.cells
      .attr("cx", d => xScale(parseFloat(d.x)))
      .attr("cy", d => yScale(parseFloat(d.y)))

    this.cells
      .attr("r", r);

    return [xScale, yScale];
  }

  // Constructing the svg
  async render(w, h) {
    const data = this.data;

    const min = getMin(w, h);

    const r = 0.003 * min;

    const _this = this;

    this.svg
      .attr("width", w)
      .attr("height", h)

    //To avoid multiples of the same shape(cell) when re-rendering
    this.gCells.selectAll("*").remove();

    //Adding the cells(circle shapes)
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

    //Tooltip functionality
    this.cells
      .on("mouseover", function (m, d) { 
        d3.select(this)
          .transition()
          .duration(100)
          .attr('r', r * 1.5) //getting bigger on hover
          .style("stroke", "black") //getting a black outline
          .style("opacity", 1)
        //Positioning tooltip
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
      //Remove all tooltip effects
      .on("mouseout", function () {
        tooltip.style("visibility", "hidden");
        d3.select(this)
          .transition()
          .duration(100)
          .attr('r', r)
          .style("stroke", "none")
          .style("opacity", 0.8)
      })


    // Pan and mouse zoom
    this.svg
      .attr('class', 'mouse-capture')
      .attr('x', -w)
      .attr('y', -h)
      .attr('width', w)
      .attr('height', h)
      .style('fill', 'white')
      .lower() //creates a layer underneath the diagram, so the mouse zoom and path work with React
      .call(zoomM)

    this.cells
      .attr("cx", d => xScale(parseFloat(d.x)))
      .attr("cy", d => yScale(parseFloat(d.y)))

  }

}
