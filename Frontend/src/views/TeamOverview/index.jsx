import React, { useState, useEffect } from 'react';
import {
  Switch, Route, useRouteMatch,
} from 'react-router-dom';
import HeaderView from 'components/general/HeaderView';
import TeamCard from 'components/teams/overview/TeamCard';
import TeamPage from 'views/TeamPage';
import styles from './teamOverview.module.css';
import queryMyTeams from 'shared/services/mock/teams';

function TeamOverview({ sidebarShown }) {
  const [teams, setTeams] = useState([]);
  useEffect(() => {
    queryMyTeams()
      .then((newTeams) => setTeams(newTeams))
      .catch((ignored) => { console.log(ignored); });
  }, [setTeams]);

  function onLeft(team) {
    setTeams(teams.filter((i) => i.id !== team.id));
  }

  const { path } = useRouteMatch();

  return (
    <Switch>
      <Route exact path={`${path}/`}>
        <HeaderView
          sidebarShown={sidebarShown}
          title="My Teams"
        >
          <div className={styles.content}>
            {teams.length === 0 ? 'No teams.' : ''}
            {teams.map((team) => (
              <div key={team.id}>
                <TeamCard
                  team={team}
                  onLeft={(t) => onLeft(t)}
                />
                <div className={styles.cardSpacing} />
              </div>
            ))}
          </div>
        </HeaderView>
      </Route>

      <Route path={`${path}/:id`}>
        <TeamPage sidebarShown={sidebarShown} />
      </Route>
    </Switch>
  );
}

export default TeamOverview;
