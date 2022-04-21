import React, { useState, useEffect } from 'react';
import { useHistory } from 'react-router-dom';
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

  const history = useHistory();
  const navigateToTeam = () => history.push(`/sequencer/teams/${team.id}`);

  const preventBubble = (evt) => {
    evt.stopPropagation();
  };

  const {
    name, description, adminIds,
  } = team;
  return (
    <div onClick={navigateToTeam} onKeyPress={navigateToTeam} role="button" tabIndex={-1}>
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
          // eslint-disable-next-line jsx-a11y/no-static-element-interactions
          <div key="leave" onClick={preventBubble} onKeyPress={preventBubble}>
            <TeamLeaveButton team={team} onLeft={onLeft} />
          </div>,
        ]}
      />
    </div>
  );
}

export default TeamCard;
