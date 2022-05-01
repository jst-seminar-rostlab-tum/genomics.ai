import React, { useState, useEffect } from 'react';
import { useHistory } from 'react-router-dom';
import SettingsIcon from '@mui/icons-material/Settings';
import Button from '@mui/material/Button';
import ListCard from 'components/general/ListCard';
import TeamLeaveButton from 'components/teams/overview/TeamLeaveButton';
import styles from './teamCard.module.css';

import getProfile from 'shared/services/profile';

function TeamCard({ team, onLeft }) {
  const [user, setUser] = useState({});
  useEffect(() => {
    getProfile()
      .then(setUser);
  }, [setUser]);

  const history = useHistory();
  const navigateToTeam = () => history.push(`/sequencer/teams/${team.id}`);

  const preventBubble = (evt) => {
    evt.stopPropagation();
  };

  const {
    name, description, adminIds,
  } = team;
  return (
    <ListCard
      title={name}
      description={description}
      nextToTitle={(
        <span className={styles.accessRightIndicator}>
          {adminIds.indexOf(user.id) !== -1 ? 'admin' : 'member'}
        </span>
      )}
      trailing={[
        adminIds.indexOf(user.id) !== -1 ? (
          <Button
            key="settings"
            endIcon={<SettingsIcon />}
            variant="outlined"
            sx={{ marginRight: '6px' }}
            onClick={navigateToTeam}
          >
            Settings
          </Button>
        ) : <div key="nothing" />,
        <div key="leave" onClick={preventBubble} onKeyPress={preventBubble} role="button" tabIndex={0}>
          <TeamLeaveButton team={team} onLeft={onLeft} />
        </div>,
      ]}
      onClick={navigateToTeam}
    />
  );
}

export default TeamCard;
