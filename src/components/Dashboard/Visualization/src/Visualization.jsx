import React from 'react';
import umapData from './data/retinal_bipolar_viz/umap_data';
import metaData from './data/retinal_bipolar_viz/meta_data';
import { UmapVisualization } from './container';
import { csv, select } from 'd3';
import './Visualization.css';

class Visualization extends React.Component {
  constructor(props) {
    super(props);
    this.myRef = React.createRef();
  }

  componentDidMount() {
    Promise.all([csv(metaData), csv(umapData)]).then(data => {
      const container = select(this.myRef.current);
      const viz = new UmapVisualization(container, data);
    });
  }

  render() {
    return (<div ref={this.myRef}></div>);
  }
}

export default Visualization;