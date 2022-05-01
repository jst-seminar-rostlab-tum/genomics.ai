import React from 'react';
import {
  Switch, Route, Redirect, useRouteMatch,
} from 'react-router-dom';
import Sidebar from '../Sidebar';
import Dashboard from 'views/Dashboard';
import TeamOverview from 'views/TeamOverview';
import InstitutionOverview from 'views/InstitutionOverview';
import UserProfile from 'views/UserProfile';
import Documentation from 'views/Documentation';
import Settings from 'views/Settings';
import Help from 'views/Help';
import styles from './dashboardContent.module.css';

const DashboardContent = () => {
  const { path, url } = useRouteMatch();

  return (
    <div>
      <Sidebar />
      <Switch>
        <Route exact path={`${path}/`}>
          <Redirect to={`${url}/dashboard`} />
        </Route>
        <Route path={`${path}/dashboard`}>
          <Dashboard />
        </Route>

        <Route path={`${path}/teams`}>
          <TeamOverview />
        </Route>

        <Route path={`${path}/institutions`}>
          <InstitutionOverview />
        </Route>

        <Route path={`${path}/users`}>
          <UserProfile />
        </Route>

        <Route path={`${path}/documentation`}>
          <Documentation />
        </Route>

        <Route path={`${path}/help`}>
          <Help />
        </Route>

        <Route path={`${path}/settings`}>
          <Settings
            className={styles.subpage}
          />
        </Route>
      </Switch>
    </div>
  );
};

export default DashboardContent;
