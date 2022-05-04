import React, { useState, useEffect } from 'react';
import {
  Switch, Route, useRouteMatch,
} from 'react-router-dom';
import Button from '@mui/material/Button';
import HeaderView from 'components/general/HeaderView';
import TeamCard from 'components/teams/TeamCard';
import TeamPage from 'views/TeamPage';
import styles from './teamOverview.module.css';
import TeamService from 'shared/services/Team.service';
import TeamCreationDialog from 'components/teams/overview/TeamCreationDialog';

function TeamOverview() {
  const [teams, setTeams] = useState([]);
  useEffect(() => {
    TeamService.getMyTeams()
      .then(setTeams)
      .catch(console.log);
  }, [setTeams]);

  function onLeft(team) {
    setTeams(teams.filter((i) => i.id !== team.id));
  }

  const [createOpen, setCreateOpen] = useState(false);

  const { path } = useRouteMatch();

  return (
    <Switch>
      <Route exact path={`${path}/`}>
        <>
          <HeaderView
            title="My Teams"
            replaceHeaderRight={(
              <Button onClick={() => setCreateOpen(true)}>Create</Button>
            )}
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
          <TeamCreationDialog
            open={createOpen}
            handleClose={() => setCreateOpen(false)}
            onCreated={(newTeam) => setTeams([...teams, newTeam])}
          />
        </>
      </Route>

      <Route path={`${path}/:id`}>
        <TeamPage />
      </Route>
    </Switch>
  );
}

export default TeamOverview;
