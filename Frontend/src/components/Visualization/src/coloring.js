import * as cons from "./constants";
import * as d3 from "d3";

export const listColoringDomain = (data, mode) => {
    let coloringDomain = data.map(x => x[mode]).filter((x, i, a) => a.indexOf(x) == i);
    return coloringDomain;
  }

export const setColoring = (mode, data) => {
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

  //Get the possible color modes with their values
export const getColoringModes = (data) =>
Object.assign({},
  ...(Object.keys(data[0])
    .map(d => d.trim())
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
      a[d] = Object.assign({},
        ...colorDomain)
      return a;
    })));