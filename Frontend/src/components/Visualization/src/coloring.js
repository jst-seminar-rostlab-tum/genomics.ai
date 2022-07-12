import * as cons from "./constants";
import * as d3 from "d3";

//Gets a list of the unique values of a given cathegory, e.g. cell_type
export const listColoringDomain = (data, mode) => {
    let coloringDomain = data.map(x => x[mode]).filter((x, i, a) => {
      return a.indexOf(x) == i 
    });
    return coloringDomain;
  }

//Makes a coloring scale, values are sorted alphabetically
export const setColoring = (mode, data) => {
    let colorScale;
    if (!parseFloat(data[0][mode])) {
      //For non-numeric values we have separate colours
      let colorDomain = listColoringDomain(data, mode).sort();
      colorScale =
        d3.scaleOrdinal()
          .domain(colorDomain)
          .range(cons.colors.slice(0, colorDomain.length));
  
    }
    else {
      //For numeric values, we have a gradient
      var colorDomain = data.map(d => parseFloat(d[mode]));
      colorScale =
        d3.scaleLinear()
          .domain([d3.min(colorDomain), d3.max(colorDomain)])
          .range(cons.gradientColors);
  
    }
    return colorScale;
  }

let keyExists = (key, obj) => key in obj;


//This function returns the coloring modes as an object in the form
//"category":
//    "subcat1": "color1"
//    "subcat2": "color2"
export const getColoringModes = (data) =>
Object.assign({},
  //Getting all the categories a cell has
  ...(Object.keys(data[0])
    .map(d => d.trim())
    //Filtering out the coordinates and categories with no name
    .filter(d => (d !== "x" && d !== "y" && d !== ""))
    .map(d => {
      var a = {};
      //If the category has numeric values
      //values and colors are assigned after sorting to match the color scale
      if (parseFloat(data[0][d])) {
        //Choose min and max
        const max = d3.max(data.map(val => parseFloat(val[d])));
        const min = d3.min(data.map(val => parseFloat(val[d])));
        //Make an object of type: <category value>: <color> 
        a[d] = Object.assign({}, ...listColoringDomain(data, d).map(val => parseFloat(val)).filter(val => val == max || val == min).sort((a, b) => a - b)
          .map((d, i) => {
            const obj = new Object();
            obj[d] = cons.gradientColors[i];
            return obj;
          }));
        return a;
      }
      //If the category has non-numeric values
      const colorDomain = listColoringDomain(data, d).sort()
      //Make an object of type: <category value>: <color> 
        .map((d, i) => {
          const obj = new Object();
          obj[d] = cons.colors[i];
          return obj;
        });
      a[d] = Object.assign({},
        ...colorDomain)
      return a;
    })));


//Unused due to non-general accuracy, aims to improve performance
// get the possible color modes with their values
export const getColoringModes2 = (data) => { // convert K-V -> key->{..}
  var colorModes = new Object();
  for(const block of data) {
    for(const [key, value] of Object.entries(block)) {
      key.trim();
      if(key != "x" && key != "y" && key != "") {
        if(!keyExists(key, colorModes) || !keyExists(value, colorModes[key])) {
          if(!keyExists(key, colorModes)) {
              colorModes[key] = {};
          }
          if(!isNaN(value)) { // expand min/max at the same time; without iterating just to figure out the min/max
            if (value < Object.keys(colorModes[key])[0] || Object.keys(colorModes[key]).length == 0) {
              delete colorModes[key][Object.keys(colorModes[key])[0]];
              colorModes[key][value] = cons.gradientColors[0];
            } else if (Number(value) > Number(Object.keys(colorModes[key])[1]) || Object.keys(colorModes[key]).length <= 1) {
              delete colorModes[key][Object.keys(colorModes[key])[1]];
              colorModes[key][value] = cons.gradientColors[1];
            }
          }else{
            colorModes[key][value] = "";
          }
        }
      }
    }
  }
  const cloneColorModes = { ...colorModes };
  for(const [key, value] of Object.entries(cloneColorModes)) {
    if(! keyExists('0', value)) {
      colorModes[key] = {};
      const tmp = Object.keys(value).sort(function (a, b) {
        return ('' + a).localeCompare(b); // force into string, not too important since we have the parseFloat but it covers edge cases
      })
      for(var i = 0; i < tmp.length; i++) {
        colorModes[key][tmp[i]] = cons.colors[i];
      }
    }
  }
  return colorModes;
}
