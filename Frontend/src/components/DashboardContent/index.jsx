import React, { useState } from 'react';
import {
  Switch, Route, Redirect, useRouteMatch,
} from 'react-router-dom';
import Sidebar from '../Sidebar';
import CommunityOverview from 'views/CommunityOverview';
import GeneMapper from 'views/GeneMapper';
import UserProfile from 'views/UserProfile';
import Documentation from 'views/Documentation';
import Settings from 'views/Settings';
import Help from 'views/Help';
import styles from './dashboardContent.module.css';
import SearchPage from 'views/SearchPage';
import { useAuth } from 'shared/context/authContext';
import Explore from 'views/Explore';
import InstitutionPage from 'views/InstitutionPage';
import TeamPage from 'views/TeamPage';
import { InstitutionProvider } from 'shared/context/institutionContext';

const DashboardContent = () => {
  const [sidebarShown] = useState(true);
  const [, setUser] = useAuth();

  const { path, url } = useRouteMatch();

  return (
    <div>
      <Sidebar setUser={setUser} />
      <Switch>
        <Route exact path={`${path}/`}>
          <Redirect to={`${url}/genemapper`} />
        </Route>

        <Route path={`${path}/genemapper`}>
          <GeneMapper sidebarShown={sidebarShown} />
        </Route>

        <Route path={`${path}/community`}>
          <InstitutionProvider>
            <CommunityOverview />
          </InstitutionProvider>
        </Route>
        <Route path={`${path}/institutions/:id`}>
          <InstitutionPage />
        </Route>
        <Route path={`${path}/teams/:id`}>
          <TeamPage />
        </Route>

        <Route path={`${path}/users/:id`}>
          <UserProfile />
        </Route>

        <Route path={`${path}/documentation`}>
          <Documentation />
        </Route>

        <Route path={`${path}/help`}>
          <Help />
        </Route>

        <Route path={`${path}/search/:searchCategory`}>
          <SearchPage sidebarShown={sidebarShown} />
        </Route>

        <Route path={`${path}/explore`} render={() => <Explore />} />

        <Route path={`${path}/settings`}>
          <Settings className={styles.subpage} />
        </Route>
      </Switch>
      {/* <Footer /> */}
    </div>
  );
};

export default DashboardContent;
