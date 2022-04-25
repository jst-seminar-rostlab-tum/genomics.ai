import React, { useState } from 'react';
import {
  Switch, Route, Redirect, useRouteMatch,
} from 'react-router-dom';
import Sidebar from '../Sidebar';
import Dashboard from 'views/Dashboard';
import TeamOverview from 'views/TeamOverview';
import InstitutionOverview from 'views/InstitutionOverview';
import Documentation from 'views/Documentation';
import Settings from 'views/Settings';
import Help from 'views/Help';
import styles from './dashboardContent.module.css';

const DashboardContent = (props) => {
  const [sidebarShown, setSidebarShown] = useState(true);
  const toggleSidebar = () => setSidebarShown(!sidebarShown);
  const { user, setUser } = props;

  const { path, url } = useRouteMatch();

  return (
    <div>
      <Sidebar
        toggleSidebar={toggleSidebar}
        sidebarShown={sidebarShown}
        user={user}
        setUser={setUser}
      />
      <Switch>
        <Route exact path={`${path}/`}>
          <Redirect to={`${url}/dashboard`} />
        </Route>
        <Route path={`${path}/dashboard`}>
          <Dashboard sidebarShown={sidebarShown} />
        </Route>

        <Route path={`${path}/teams`}>
          <TeamOverview sidebarShown={sidebarShown} />
        </Route>

        <Route path={`${path}/institutions`}>
          <InstitutionOverview sidebarShown={sidebarShown} />
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
      </Switch>
    </div>
  );
};

export default DashboardContent;
