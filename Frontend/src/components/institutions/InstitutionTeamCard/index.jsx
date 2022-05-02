import React from 'react';
import { useHistory } from 'react-router-dom';
import ListCard from 'components/general/ListCard';
import styles from './institutionTeamCard.module.css';

function InstitutionTeamCard({ team /* onLeft, institution */ }) {
  const history = useHistory();
  const navigateToTeam = () => history.push(`/sequencer/teams/${team.id}`);

  const {
    name, description, visibility,
  } = team;

  return (
    <ListCard
      title={name}
      description={description}
      nextToTitle={(
        <span className={styles.accessRightIndicator}>
          {visibility}
        </span>
      )}
      onClick={navigateToTeam}
    />
  );
}

export default InstitutionTeamCard;
