import React, { useState, useEffect } from 'react';
import { useHistory } from 'react-router-dom';
import SettingsIcon from '@mui/icons-material/Settings';
import Button from '@mui/material/Button';
import ListCard from 'components/general/ListCard';
import TeamLeaveButton from 'components/teams/overview/TeamLeaveButton';
import styles from './institutionTeamCard.module.css';

import getUser from 'shared/services/mock/user';

function InstitutionTeamCard({ team, onLeft }) {
  const [user, setUser] = useState({});
  useEffect(() => {
    getUser()
      .then(setUser);
  }, [setUser]);

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
