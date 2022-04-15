import React from 'react';
import {
  Switch, Route, Redirect, useRouteMatch,
} from 'react-router-dom';
import Dashboard from '../Dashboard';
import Documentation from '../Pages/Documentation/Documentation';
import Settings from '../Pages/Settings/Settings';
import Help from '../Pages/Help/Help';
import styles from './dashboardContent.module.css';
import Sidebar from '../Sidebar/Sidebar';

const DashboardContent = (props) => {
  const { user, setUser } = props;

  const { path, url } = useRouteMatch();

  return (
    <>
      <Sidebar setUser={setUser} />
      <Switch>
        <Route exact path={`${path}/`}>
          <Redirect to={`${url}/dashboard`} />
        </Route>
        <Route path={`${path}/dashboard`}>
          <Dashboard />
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
            user={user}
            setUser={setUser}
          />
        </Route>
      </Switch>
    </>
  );
};

export default DashboardContent;
