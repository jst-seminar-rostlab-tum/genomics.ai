import * as d3 from "d3";
import * as cons from "./constants";
import {initZoom} from "./newZoom"
//TODO: Fix constants
//TODO: Convert x and y to float beforehand before to_csv
//TODO: Add a ref attribute

const originalOpacity = 0.75;

const listColoringDomain = (data, mode) => {
    let coloringDomain = data.map(x => x[mode]).filter((x, i, a) => a.indexOf(x) == i);
    return coloringDomain;
 }


 //For Before and After query
 const queryAfter = (cells) =>{
    cells
    .style("opacity", d => {
      if (!d.ref){
        return 0;
      }
      else{
        return originalOpacity;
      }
    })
  
  }
  
  const queryBefore = (smth) =>{
    cells
    .style("opacity", originalOpacity)
  }

  const setColoring = (mode, data) => {
    if(mode === "cell_type" || mode === "batch")
    {
      let colorDomain = listColoringDomain(data, mode);
      var colorScale = 
      d3.scaleOrdinal()
      .domain(...colorDomain)
      .range(cons.colors.slice(0,colorDomain.length));
    }
    else{
      var colorDomain = data.map(d => parseFloat(d[mode]));
      var colorScale = 
      d3.scaleLinear()
      .domain([d3.min(colorDomain),d3.max(colorDomain)])
      .range(["#42e91f", "#0a43c8"]);
    
    }
    return colorScale;
  }

  const addGroup = (svg, id) => {
    return svg.append('g').attr('id', `${id}`);
  };

  export class UmapVisualization2 {


    constructor(container, data) {
  
      this.svg = d3.select(container).append('svg');
      this.gCells = addGroup(this.svg, 'cells');
      this.gLabels = addGroup(this.svg, 'labels');
      this.data = data;
    };

    setColorMode(mode) {
      const colorScale = setColoring(mode, this.data);
      this.cells = this.cells.style("fill", (d)=> colorScale(d[mode]));
    }

    async render(w, h){
        const data = this.data; //dataUnpacked;
        
        //Scales
        const  xScale = d3.scaleLinear()
        .domain([d3.min(data.map(d => parseFloat(d.x))), d3.max(data.map(d => parseFloat(d.x)))])
        .range([cons.margin, w-cons.margin]);
  
        const  yScale = d3.scaleLinear()
        .domain([d3.min(data.map(d => parseFloat(d.y))), d3.max(data.map(d => parseFloat(d.y)))])
        .range([cons.margin, h-cons.margin]);
 
         //svg
         this.svg
         .attr("width", w)
         .attr("height", h)
         
 
         //Circle cells
         this.cells = this.gCells
         .selectAll("cell")
         .data(data)
         .enter()
         .append("circle")
         .attr("cx", d => xScale(parseFloat(d.x)))
         .attr("cy", d => yScale(parseFloat(d.y)))
         .attr("r", cons.initialPointRadius)
        //  .style("fill", (d)=> colorScale(d[mode]))
         .style("fill", "black")
         .style("opacity", originalOpacity)

         //Pan and mouse zoom
         //initZoom();

     }

    }
 