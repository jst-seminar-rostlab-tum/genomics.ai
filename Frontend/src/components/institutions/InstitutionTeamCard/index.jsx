import React from 'react';
import { useHistory } from 'react-router-dom';
import ListCard from 'components/general/ListCard';
import TeamLeaveButton from 'components/teams/overview/TeamLeaveButton';
import styles from './institutionTeamCard.module.css';
import { useAuth } from 'shared/context/authContext';

function InstitutionTeamCard({ team, onLeft }) {
  const [user] = useAuth();

  const history = useHistory();
  const navigateToTeam = () => history.push(`/sequencer/teams/${team.id}`);

  const preventBubble = (evt) => {
    evt.stopPropagation();
  };

  const {
    name, description, adminIds, visibility
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
      trailing={[
        adminIds.indexOf(user.id) !== -1 ? (
          <div key="leave" onClick={preventBubble} onKeyPress={preventBubble} role="button" tabIndex={0}>
            <TeamLeaveButton team={team} onLeft={onLeft} />
          </div>
        ) : <div key="nothing" />,
      ]}
      onClick={navigateToTeam}
    />
  );
}

export default InstitutionTeamCard;
