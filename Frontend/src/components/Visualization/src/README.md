# Visualization
The visualization is based on the class called UmapVisualization2 and all of the key elements for the end product are contained in there.
The files used for the GeneMapper are: umapVisualization.js, coloring.js, newZoom.js, barChart.js for the functionality and constants.js, bar.css, tooltip.css for the styling. The others belong to the GeneCruncher.

***
## UmapVisualization2
This class has the following constructor:
```
constructor(container, data, containerBar) {
    d3.select(container).selectAll("*").remove();
    this.svg = d3.select(container).append('svg');
    this.gCells = addGroup(this.svg, 'cells');
    this.coloringModes = getColoringModes(data);
    this.tooltip = d3.select(container).append("div");
    this.mode = undefined;
    this.barChartBatch = addBarPlotBatch(containerBar, data);
    if(Object.keys(this.coloringModes).includes("cell_type")){
      d3.select(containerBar).append("div");
      this.barChartCell = addBarPlotCell(containerBar, data);
    }
    this.data = data;
  };
```
### Constructor:
**Arguments:**
* _container_ - the container in which the svg(scatter plot) is. On the frontend it is a div
* _data_ - is the already parsed .csv data
* _containerBar_ - the container in which the additional diagrams can be found. Again a div on the frontend

**Class Attributes**
* _gCells_ - cell group
* _cells_ - the cells inside the gCells group
* _coloringModes_ - every category with its sub-categories in the form \{category1: \{subCategory1: color\} \.\.\. \}, sorted alphabetically
* _mode_ - can take a value from the set of attributes
* _colorScale_ - the color scale that provides the plot coloring according to the mode, initial color, however, black

### Class functions
* _before_ and _after_ - responsible for the on and off toggle for the sub-categories
* _setColorMode_ - changes the colorScale according to the color mode chosen
* _resize_ - sets the dimension scales for the plot, the dot coordinates and the dot radius
*_render_ - visualized the data and provides the interactive tooltip

***
## newZoom
Provides all zoom functionalities
* _zoomM_ - the scale for the zoom functionalities and responsible for mouse zoom and pat
* _zoomInN_, _zoomOutN_, _resetZoom_ - resposible for zooming in, out and resetting the zoom scale

***
## coloring
Provides all coloring functionalities
* _listColoringDomain_ - returns a list of the sub-categories
* _setColoring_ - produces the colorScale for _setColorMode_ 

***
## barChart
Provides the functionalities for all additional diagrams
* _groupBy_ - groups the data according to the given _cat_, category, attribute
