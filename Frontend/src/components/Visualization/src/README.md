# Visualization
Generally is the visualization responsible for the logic and drawing (not the placement) of the umap and the bar graphs. The visualizing is done with d3.
Most elements are done with possible extensions in mind.
The main class is called UmapVisualization2 and all of the key elements for the end product are contained in there.
The files used for the GeneMapper are: umapVisualization.js, coloring.js, newZoom.js, barChart.js for the functionality and constants.js, bar.css, tooltip.css for the styling. The others belong to the GeneCruncher.

***
## UmapVisualization2
This class has the following constructor:
```
class UmapVisualization2 {

  constructor(container, data) {
    d3.select(container).selectAll("*").remove();
    this.svg = d3.select(container).append('svg')
      .attr('id', 'vis_svg');
    this.gCells = addGroup(this.svg, 'cells');
    this.coloringModes = getColoringModes(data);
    this.tooltip = d3.select(container).append("div");
    this.mode = undefined;
    this.data = data;
    this.graphs = ["batch", "cell"];
    this.hiddenCells = [];
  };
```
### Constructor:
**Arguments:**
* _container_ - the container in which the svg(scatter plot) is. On the frontend it is a div
* _data_ - is the already parsed .csv data

**Class Attributes**
* _gCells_ - cell group
* _cells_ - the cells inside the gCells group
* _coloringModes_ - every category with its sub-categories in the form \{category1: \{subCategory1: color\} \.\.\. \}, sorted alphabetically
* _mode_ - can take a value from the set of attributes
* _colorScale_ - the color scale that provides the plot coloring according to the mode, initial color, however, black
* _graphs_ - the available graphs (can be extracted via getAvailableGraphs())
* _hiddenCells_ - contains a list of [category, value] which specifies the currently hidden cells

### Class functions
* _drawGraph_ - draws a specified graph with a given height and width in a given container (only graphs from _graphs_ are possible)
* _before_ and _after_ - responsible for the on and off toggle for the sub-categories
* _beforeAll_ - resets the visibility of all cells
* _showReference_, _hideReference_, _showQuery_ and _hideQuery_  - are extra methods that use the _before_/_after_ methods
* _reduceOpacity_ - reduces the opacity of a category and value (only one at the same time is supported at the moment)
* _resetOpacity_ - sets the opacity of all cells to the original value
* _predictedCellsTransparten_ and _predictedCellsVisible_ - use _reduceOpacity_ and _resetOpacity_ to provide the functionality for the predicted cells button
* _setColorMode_ - changes the colorScale according to the color mode chosen
* _resize_ - sets the dimension scales for the plot, the dot coordinates and the dot radius
*_render_ - visualized the data and provides the interactive tooltip

***
## newZoom
Provides all zoom functionalities
* _zoomM_ - the scale for the zoom functionalities and responsible for mouse zoom and pat
* _zoomInN_, _zoomOutN_, _resetZoom_ - resposible for zooming in, out and resetting the zoom scale  
-> the mouse wheel zoom is called in umapVisualization, the button zoom is called in the various front end files  
-> is called on the specific id of the svg, because it can come to problems otherwise  

***
## coloring
Provides all coloring functionalities
* _listColoringDomain_ - returns a list of the sub-categories
* _setColoring_ - produces the colorScale for _setColorMode_ 

***
## barChart
Provides the functionalities for all additional diagrams
* _groupBy_ - groups the data according to the given _cat_, category, attribute
* _addBarPlot_ - draws the chart in the given barContainer, with the other given attributes
