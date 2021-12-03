import React, { useState, useEffect } from 'react';
import {
  Switch, Route, Redirect, useRouteMatch,
} from 'react-router-dom';
import Sidebar from '../Sidebar/Sidebar';
import Dashboard from '../Dashboard';
import NavigationBar from '../NavigationBar/NavigationBar';
import Documentation from '../Pages/Documentation/Documentation';
import Settings from '../Pages/Settings/Settings';
import Help from '../Pages/Help/Help';
import styles from './dashboardContent.module.css';
import VisualizationPage from '../Pages/VisualizationPage/VisualizationPage';

const DashboardContent = (props) => {
  const [sidebarShown, setSidebarShown] = useState(true);
  const toggleSidebar = () => setSidebarShown(!sidebarShown);
  const { user, setUser } = props;

  const { path, url } = useRouteMatch();

  // TODO: fetch the visualization results here, currently setting the results to a static array: [1]
  const visualizationResults = [1];

  return (
    <Switch>
      <div>
        <NavigationBar sidebarShown={sidebarShown} user={user} setUser={setUser} />
        <Sidebar
          toggleSidebar={toggleSidebar}
          sidebarShown={sidebarShown}
        />
        <Route exact path={`${path}/`}>
          <Redirect to={`${url}/dashboard`} />
        </Route>
        <Route path={`${path}/dashboard`}>
          <Dashboard sidebarShown={sidebarShown} />
        </Route>

        <Route path={`${path}/documentation`}>
          <Documentation sidebarShown={sidebarShown} />
        </Route>

        <Route path={`${path}/help`}>
          <Help sidebarShown={sidebarShown} />
        </Route>

        <Route path={`${path}/settings`}>
          <Settings
            className={sidebarShown ? styles.subpage : styles.subpageSidebarCollapsed}
            sidebarShown={sidebarShown}
          />
        </Route>
         {/* TODO: Looping over all visualization projects
          Check backend specification to determine how to create the route for it
        TODO: add the id to the visualization page
        */}
        {
          visualizationResults.map((id) => (
            <Route path={`${path}/result${id}`}>
              <VisualizationPage id={id} />
            </Route>
          ))
        }
        {/*The lines below are for testing*/}
        {/* <Route path={`${path}/visualization`}>
          <VisualizationPage id={id} />
        </Route> */}
      </div>
    </Switch>
  );
};

export default DashboardContent;
