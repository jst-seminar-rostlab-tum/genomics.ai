import React, { useState, useEffect } from 'react';
import HeaderView from 'components/general/HeaderView';
import TeamCard from 'components/teamOverview/TeamCard';
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

  return (
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
  );
}

export default TeamOverview;
