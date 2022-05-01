import React from 'react';
import { useHistory } from 'react-router-dom';
import ListCard from 'components/general/ListCard';
import TeamLeaveButton from 'components/teams/overview/TeamLeaveButton';
import styles from './institutionTeamCard.module.css';
import { useAuth } from 'shared/context/authContext';

function InstitutionTeamCard({ team, onLeft, institution }) {
  const [user] = useAuth();

  const history = useHistory();
  const navigateToTeam = () => history.push(`/sequencer/teams/${team.id}`);

  const preventBubble = (evt) => {
    evt.stopPropagation();
  };

  function isAdmin() {
    return (institution.adminIds || []).includes(user._id);
  }

  const {
    name, description, visibility
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
