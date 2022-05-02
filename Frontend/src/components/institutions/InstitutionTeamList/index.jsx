import React, { useState, useEffect } from 'react';
import CircularProgress from '@mui/material/CircularProgress';
import styles from './institutionTeamList.module.css';
import TeamService from 'shared/services/Team.service';
import InstitutionTeamCard from 'components/institutions/InstitutionTeamCard';

function InstitutionTeamList({ onLeft, institution }) {
  const [teams, setTeams] = useState([]);
  const [teamsLoaded, setTeamsLoaded] = useState(false);

  useEffect(async () => {
    setTeams(await TeamService.getInstitutionTeams(institution.id));
    setTeamsLoaded(true);
  }, []);

  if (!teamsLoaded) {
    return <CircularProgress />;
  }

  return (
    <div className={styles.content}>
      {teams.length === 0 ? 'No teams.' : ''}
      {teams.map((team) => (
        <div key={team.id}>
          <InstitutionTeamCard
            team={team}
            onLeft={onLeft}
            institution={institution}
          />
          <div className={styles.cardSpacing} />
        </div>
      ))}
    </div>
  );
}

export default InstitutionTeamList;