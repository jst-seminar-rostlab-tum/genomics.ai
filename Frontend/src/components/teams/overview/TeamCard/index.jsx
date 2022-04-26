import React, { useState, useEffect } from 'react';
import SettingsIcon from '@mui/icons-material/Settings';
import Button from '@mui/material/Button';
import ListCard from 'components/general/ListCard';
import TeamLeaveButton from 'components/teams/overview/TeamLeaveButton';
import styles from './teamCard.module.css';

import getUser from 'shared/services/mock/user';

function TeamCard({ team, onLeft }) {
  const [user, setUser] = useState({});
  useEffect(() => {
    getUser()
      .then(setUser);
  }, [setUser]);

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
          >
            Settings
          </Button>
        ) : <div key="nothing" />,
        <TeamLeaveButton key="leave" team={team} onLeft={onLeft} />,
      ]}
    />
  );
}

export default TeamCard;
