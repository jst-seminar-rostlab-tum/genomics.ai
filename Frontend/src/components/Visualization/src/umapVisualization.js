import * as d3 from "d3";
import * as cons from "./constants";
import { zoomM } from "./newZoom";
import "./tooltip.css";

//TODO: Fix constants
//TODO: Refactor coloring functions

//
//Bugs:
//AttrAdding doesn't work properly
//Cells are overwritten - *4
//

const getMin = (w, h) => {
  return d3.min([h, w]);
}

//Hide (for every possible attribute)
const queryAfter = (cells, category, value) => {
  cells
  .style("visibility", function(d){
    // console.log((d3.select(this).attr("active")));
    console.log("here");
    if (d[category] == value){
      d3.select(this)
      .attr("active", false);
      return "hidden";
    }
  });

}

//Show
const queryBefore = (cells) => {
  cells
  .style("visibility", (d) => (d[category]===value? "hidden" : "visible"))
}

//Reference before and after


const listColoringDomain = (data, mode) => {
  let coloringDomain = data.map(x => x[mode]).filter((x, i, a) => a.indexOf(x) == i);
  return coloringDomain;
}

//Create a color scale based on the chosen mode and its values
const setColoring = (mode, data) => {
  let colorScale;
  if (!parseFloat(data[0][mode])) {
    let colorDomain = listColoringDomain(data, mode).sort();
    colorScale =
      d3.scaleOrdinal()
        .domain(colorDomain)
        .range(cons.colors.slice(0, colorDomain.length));

  }
  else {
    var colorDomain = data.map(d => parseFloat(d[mode]));
    colorScale =
      d3.scaleLinear()
        .domain([d3.min(colorDomain), d3.max(colorDomain)])
        .range(cons.gradientColors);

  }
  return colorScale;
}

const addGroup = (svg, id) => {
  return svg.append('g').attr('id', `${id}`);
};

// //TODO: Fix
// const attrAdding = (cells, attrDom) => {
//   for (let i = 0; i < attrDom.length; i++) {
//     var att = attrDom[i];
//     cells
//       .attr(att, (d) => { return d[att] });
//   }

// }

//Get the possible color modes with their values
const getColoringModes = (data) =>
  Object.assign({},
    ...(Object.keys(data[0])
      .filter(d => (d !== "x" && d !== "y" && d !== ""))
      .map(d => {
        var a = {};
        if (parseFloat(data[0][d])) {
          const max = d3.max(data.map(val => parseFloat(val[d])));
          const min = d3.min(data.map(val => parseFloat(val[d])));
          a[d] = Object.assign({}, ...listColoringDomain(data, d).map(val => parseFloat(val)).filter(val => val == max || val == min).sort((a, b) => a - b)
            .map((d, i) => {
              const obj = new Object();
              obj[d] = cons.gradientColors[i];
              return obj;
            }));
          return a;
        }
        const colorDomain = listColoringDomain(data, d).sort()
          .map((d, i) => {
            const obj = new Object();
            obj[d] = cons.colors[i];
            return obj;
          });
        // console.log(colorDomain);
        a[d] = Object.assign({},
          ...colorDomain)
        return a;
      })));

export class UmapVisualization2 {


  constructor(container, data) {
    d3.select(container).selectAll("*").remove();
    this.svg = d3.select(container).append('svg');
    this.gCells = addGroup(this.svg, 'cells');
    this.gLabels = addGroup(this.svg, 'labels');
    this.coloringModes = getColoringModes(data);
    this.tooltip = d3.select(container).append("div");
    this.mode = undefined;
    this.data = data;
  };

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
    const min = getMin(w, h);
    const data = this.data;

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
      .attr("active", true)

      // queryAfter(this.cells, "cell_type", "Pancreas Beta")
    
      // queryAfter(this.cells, "cell_type", "Pancreas Alpha")

    //Addiding
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
    // .on('click', function (d, i){
    //     const showMore = document.getElementById("moreVisul");
    //     if(showMore.style.display === 'none') {
    //         showMore.style.display = 'block';
    //     } else {
    //         showMore.style.display = 'none';
    //     }
    // });


    //Pan and mouse zoom
    this.svg
      .attr('class', 'mouse-capture')
      .attr('x', -w)
      .attr('y', -h)
      .attr('width', w)
      .attr('height', h)
      .style('fill', 'white')
      .lower()
      .call(zoomM);

        console.log(this.cells)
  }

}
