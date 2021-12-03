import React, { useState } from 'react';
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
  const visualizationResults = [1];

  return (
    <div>
      <NavigationBar sidebarShown={sidebarShown} user={user} setUser={setUser} />
      <Sidebar
        toggleSidebar={toggleSidebar}
        sidebarShown={sidebarShown}
      />

      <Switch>
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
            user={user}
            setUser={setUser}
          />
        </Route>

        {/* Looping over all visualization projects
        Check backend specification to determine how to create the route for it
        */}
        {
          visualizationResults.map((id) => (
            <Route path={`/result${id}`} key={id}>
              <VisualizationPage />
            </Route>
          ))
        }
      </Switch>

    </div>

  );
};

export default DashboardContent;
