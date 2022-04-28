import { ContentCutOutlined } from "@mui/icons-material";
import * as d3 from "d3";
import { selectAll } from "d3";
import * as cons from "./constants";
import {zoomM} from "./newZoom"
//TODO: Fix constants
//TODO: Convert x and y to float beforehand before to_csv
//TODO: Add a ref attribute


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
        return cons.originalOpacity;
      }
    })
  
  }
  
  const queryBefore = (smth) =>{
    cells
    .style("opacity", cons.originalOpacity)
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
      d3.select(container).selectAll("*").remove();
      this.svg = d3.select(container).append('svg');
      this.gCells = addGroup(this.svg, 'cells');
      this.gLabels = addGroup(this.svg, 'labels');
      this.coloringModes = Object.assign({},
      ...(Object.keys(data[0])
      .filter(d => (d != "x" && d != "y" && d != ""))
      .map(d =>{ var a = {};
                  if (!parseFloat(data[0][d]))
                  {
                    a[d] = [d3.min(data.map(value => parseFloat(value[d]))), d3.max(data.map(value => parseFloat(value[d])))];
                    return a[d] ;
                  }
                  a[d] = listColoringDomain(data,d)
                return a;} )));
      this.data = data;
    };

    setColorMode(mode) {
      const colorScale = setColoring(mode, this.data);
      this.cells = this.cells.style("fill", (d)=> colorScale(d[mode]));
    }

    //Width and height should be the size of the container, not square
    async render(w, h){
      console.log(this.coloringModes);
        const data = this.data; //dataUnpacked;
        
        const min = d3.min([h,w]);
        const r = 0.003*min;

        const  xScaleHelper = d3.scaleLinear()
        .domain([d3.min(data.map(d => parseFloat(d.x))), d3.max(data.map(d => parseFloat(d.x)))])
        .range([cons.margin, min-cons.margin]);

        const translate = (w - cons.margin - (xScaleHelper(d3.max(data.map(d => parseFloat(d.x))))))/2;

        //Scales
        const  xScale = d3.scaleLinear()
        .domain([d3.min(data.map(d => parseFloat(d.x))), d3.max(data.map(d => parseFloat(d.x)))])
        .range([translate, min-cons.margin + translate]);
  
        const  yScale = d3.scaleLinear()
        .domain([d3.min(data.map(d => parseFloat(d.y))), d3.max(data.map(d => parseFloat(d.y)))])
        .range([cons.margin, min-cons.margin]);
 
         //svg
         this.svg
         .attr("width", w)
         .attr("height", h)
         
         this.cells = this.gCells.selectAll("*").remove();
         //Circle cells
         this.cells = this.gCells
         .selectAll("cell")
         .data(data)
         .enter()
         .append("circle")
         .attr("cx", d => xScale(parseFloat(d.x)))
         .attr("cy", d => yScale(parseFloat(d.y)))
         .attr("r", r)
        //  .style("fill", (d)=> colorScale(d[mode]))
         .style("fill", cons.fill)
         .style("opacity", cons.originalOpacity)
         
         
         //Pan and mouse zoom
         
         this.svg
            .attr('class', 'mouse-capture')
            .attr('x', -w)
            .attr('y', -h)
            .attr('width', w)
            .attr('height', h)
            .style('fill', 'white')
            .lower() // put it below the map
            .call(zoomM);
          

     }

    }
 