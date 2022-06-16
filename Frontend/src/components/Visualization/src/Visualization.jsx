import React from 'react';
import { tsv, select } from 'd3';
import umapData from './data/retinal_bipolar_viz/umap_data';
import { UmapVisualization } from './container';
import { width, height } from './constants';
import './visualization.module.css';

class Visualization extends React.Component {
  constructor(props) {
    super(props);
    this.myRef = React.createRef();
    this.url = props.url;
  }

  componentDidMount() {
    tsv(this.url).then(data => {
      const container = select(this.myRef.current);
      const viz = new UmapVisualization(container, data);
      viz.render([width, height]);
    });
  }

  render() {
    return (<div ref={this.myRef}></div>);
  }
}

export default Visualization;