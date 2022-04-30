import * as d3 from "d3";
import * as cons from "./constants";
import {zoomM, zoomInN, zoomOutN} from "./newZoom"
//TODO: Fix constants
//TODO: Convert x and y to float beforehand before to_csv
//TODO: Add a ref attribute
//TODO: Refactor coloring functions


const getMin = (w,h) =>
{
  return d3.min([h,w]);
}

 //Hide (for every possible attribute)
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
  
  //Show
  const queryBefore = (cells) =>{
    cells = cells
    .style("opacity", cons.originalOpacity)
  }

  //Reference before and after
  

  const listColoringDomain = (data, mode) => {
    let coloringDomain = data.map(x => x[mode]).filter((x, i, a) => a.indexOf(x) == i);
    return coloringDomain;
 }

  //Create a color scale based on the chosen mode and its values
  const setColoring = (mode, data) => {
    if(!parseFloat(data[0][mode]))
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

  const attrAdding = (cells, attrDom) =>
  {
    for(let i=0; i < attrDom.length; i++){
      var att = attrDom[i];
      cells =
      cells
      .attr(att, (d) => d[att]);
    }
    
  }

  //Get the possible color modes with their values
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

    //Set a color mode
    setColorMode(mode) {
      const colorScale = setColoring(mode, this.data);
      if (this.cells != null){
        this.cells = this.cells.style("fill", (d)=> colorScale(d[mode]));
      }
      this.mode = mode;
    }

    //Sizing and resizing the cells
    resize(w, h) {
      const min = getMin(w,h);
      const data = this.data;

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

      this.cells = 
      this.cells
      .attr("cx", d => xScale(parseFloat(d.x)))
      .attr("cy", d => yScale(parseFloat(d.y)))


      return [xScale, yScale];
    }

    //Constructing the svg
    async render(w, h){
      
        const data = this.data; //dataUnpacked;
        
        const min = getMin(w,h);

        const r = 0.003*min;

         //svg
         this.svg
         .attr("width", w)
         .attr("height", h)
         
         //Circle cells
         this.cells = this.gCells.selectAll("*").remove();

         this.cells = this.gCells
         .selectAll("cell")
         .data(data)
         .enter()
         .append("circle")
        //  .attr("cx", d => xScale(parseFloat(d.x)))
        //  .attr("cy", d => yScale(parseFloat(d.y)))
         .attr("r", r)
         .style("fill", (d) => {
           if (this.mode != null){
             const colorScale = setColoring(this.mode, this.data);
             return colorScale(d[this.mode]);
           }
           return cons.fill;
          })
         .style("opacity", cons.originalOpacity)
         
         //Adding attributes from the different categories for each cell
         const attrDom = Object.keys(this.coloringModes)
         attrDom.push("x", "y");
         attrAdding(this.cells, attrDom);

         this.resize(w,h);
         
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
          

     }

    }
 