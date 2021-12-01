import React, { useState } from 'react';
import { BrowserRouter as Router, Switch, Route } from 'react-router-dom';
import Sidebar from '../Sidebar/Sidebar';
import Dashboard from '../Dashboard';
import NavigationBar from '../NavigationBar/NavigationBar';
import Documentation from '../Pages/Documentation/Documentation';
import Settings from '../Pages/Settings/Settings';
import Help from '../Pages/Help/Help';
import VisualizationPage from '../Pages/VisualizationPage/VisualizationPage';
import styles from './dashboardContent.module.css';

const DashboardContent = () => {
  const [sidebarShown, setSidebarShown] = useState(true);
  const toggleSidebar = () => setSidebarShown(!sidebarShown);
  const [visualizationResults, setVisualizationResults] = useState([1]);

  //concatenate the visualization result project id with the route for it

  //testing the visualization results

  return (
    <Router>
      <NavigationBar sidebarShown={sidebarShown} />
      <Sidebar
        toggleSidebar={toggleSidebar}
        sidebarShown={sidebarShown}
      />

      <Switch>
        <Route path="/dashboard">
          <Dashboard sidebarShown={sidebarShown} />
        </Route>

        <Route path="/documentation">
          <Documentation sidebarShown={sidebarShown} />
        </Route>

        <Route path="/help">
          <Help sidebarShown={sidebarShown} />
        </Route>

        <Route path="/settings">
          <Settings
            className={sidebarShown ? styles.subpage : styles.subpageSidebarCollapsed}
            sidebarShown={sidebarShown}
          />
        </Route>

        {/*Looping over all visualization projects
        Check backend specification to determine how to create the route for it
        */}
        {
          visualizationResults.map((id) => {
            return (<Route path={"/result" + id}>
              <VisualizationPage />
            </Route>)
          })
        }
      </Switch>
    </Router>
  );
};

export default DashboardContent;
