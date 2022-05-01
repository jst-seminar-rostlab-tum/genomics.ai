import React, { useState, useEffect } from 'react';
import CircularProgress from '@mui/material/CircularProgress';
import styles from './institutionMemberList.module.css';
import { getInstitutionTeams } from 'shared/services/mock/teams';
import TeamCard from 'components/teams/TeamCard';

function InstitutionTeamList({ institution }) {
  const [teams, setTeams] = useState([]);
  const [teamsLoaded, setTeamsLoaded] = useState(false);

  useEffect(() => {
    getInstitutionTeams(institution.id).then((loadedTeams) => {
      setTeams(loadedTeams);
      setTeamsLoaded(true);
    });
  });

  if (!teamsLoaded) {
    return <CircularProgress />;
  }

  return (
    <div className={styles.content}>
      {teams.length === 0 ? 'No teams.' : ''}
      {teams.map((team) => (
        <div key={team.id}>
          <TeamCard
            team={team}
            onLeft={(leftTeam) => {
              setTeams(teams.filter((t) => t.id !== leftTeam.id));
            }}
          />
          <div className={styles.cardSpacing} />
        </div>
      ))}
    </div>
  );
}

export default InstitutionTeamList;
