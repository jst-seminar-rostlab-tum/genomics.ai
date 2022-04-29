import * as d3 from "d3";
import * as cons from "./constants";
import {zoomM, zoomInN, zoomOutN} from "./newZoom"
//TODO: Fix constants
//TODO: Convert x and y to float beforehand before to_csv
//TODO: Add a ref attribute
//TODO: Refactor coloring functions

const listColoringDomain = (data, mode) => {
    let coloringDomain = data.map(x => x[mode]).filter((x, i, a) => a.indexOf(x) == i);
    return coloringDomain;
 }


 //For Before and After query
 const queryAfter = (cells, category, value) =>{
    cells
    .style("opacity", d => {
      if (d[category] == value){
        return 0;
      }
      else{
        return cons.originalOpacity;
      }
    })
  
  }
  
  const queryBefore = (cells) =>{
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
      .range(cons.gradientColors);
    
    }
    return colorScale;
  }

  const addGroup = (svg, id) => {
    return svg.append('g').attr('id', `${id}`);
  };


  const getColoringModes = (data) =>
  Object.assign({},
    ...(Object.keys(data[0])
    .filter(d => (d != "x" && d != "y" && d != ""))
    .map(d =>{ var a = {};
                if (parseFloat(data[0][d]))
                {
                  let max = d3.max(data.map(val => parseFloat(val[d])));
                  let min = d3.min(data.map(val => parseFloat(val[d])));
                  a[d] = Object.assign({}, ...listColoringDomain(data,d).map(val => parseFloat(val)).filter( val => val == max || val == min).sort((a,b) => a - b)
                  .map((d,i) => {
                  const obj = new Object();
                  obj[d] = cons.gradientColors[i];
                  return obj;}));
                  return a; 
                }
                a[d] = Object.assign({},
                ...listColoringDomain(data,d)
                .map((d,i) => {
                  const obj = new Object();
                  obj[d] = cons.colors[i];
                  return obj;
                }))
              return a;} )));

  export class UmapVisualization2 {


    constructor(container, data) {
      d3.select(container).selectAll("*").remove();
      this.svg = d3.select(container).append('svg');
      this.label = this.svg.append("label");
      this.gCells = addGroup(this.svg, 'cells');
      this.gLabels = addGroup(this.svg, 'labels');
      this.coloringModes = getColoringModes(data);
      this.mode = null;
      this.data = data;
    };

    setColorMode(mode) {
      const colorScale = setColoring(mode, this.data);
      this.cells = this.cells.style("fill", (d)=> colorScale(d[mode]));
      this.mode = mode;
    }

    //Width and height should be the size of the container, not square
    async render(w, h){
      
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
         .style("fill", (d) => {
           if (this.mode != null){
             const colorScale = setColoring(this.mode, this.data);
             return colorScale(d[this.mode]);
           }
           return cons.fill;
          })
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
 